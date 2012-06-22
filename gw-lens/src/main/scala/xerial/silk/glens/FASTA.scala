//--------------------------------------
//
// FASTA.scala
// Since: 2012/06/22 10:00 AM
//
//--------------------------------------

package xerial.silk.glens

import xerial.silk.util.io.{BufferedScanner, FileSource}
import java.io._
import scala.Some
import annotation.tailrec
import xerial.silk.util.Logger
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import java.util.zip.GZIPInputStream


/**
 * FASTA file reader
 *
 * @author leo
 */
object FASTA extends Logger {

  def extractSequenceNameFrom(descriptionLine: String) = {
    val trimmed = descriptionLine.substring(1).dropWhile(c => Character.isWhitespace(c))
    trimmed.takeWhile(c => !Character.isWhitespace(c) && c != '|')
  }

  /**
   * Reading
   * @param in
   */
  private[glens] class SequenceStream(in: LineReader) {
    private var finishedReading = false
    private var initialized = false

    private def noMoreData: Boolean = {
      if (finishedReading)
        true
      else {
        val ch = in.LA1
        if (ch == '>' || ch == BufferedScanner.EOF)
          finishedReading = true
        finishedReading
      }
    }

    def isValidStream = !finishedReading

    private[FASTA] def invalidate {
      while (!noMoreData)
        in.nextLine
    }

    def toStream: Stream[String] = {
      if (initialized)
        sys.error("Cannot create stream more than once")

      initialized = true

      if (finishedReading)
        sys.error("Cannot create a stream because this stream was invalidated when reading the following FASTASequences")

      def loop: Stream[String] =
        if (noMoreData)
          Stream.empty[String]
        else
          in.nextLine.toString #:: loop
      loop
    }
  }


  trait LineReader extends Closeable {
    def LA1: Int

    def nextLine: CharSequence

    def lineCount: Int

    def close: Unit
  }


  /**
   * Add a line counter to the BufferedScanner
   * @param in
   */
  private class StandardLineReader(in: BufferedScanner) extends LineReader {
    private var _lineCount = 0

    def LA1 = in.LA(1)

    def nextLine: CharSequence = {
      _lineCount += 1
      in.nextLine
    }

    def lineCount = _lineCount

    def close = in.close
  }

  val DEFAULT_BUFFER_SIZE = 4 * 1024 * 1024 // 4MB

  private class TarFileLineReader(in: InputStream, bufferSize:Int = DEFAULT_BUFFER_SIZE) extends LineReader {
    private val tarIn = new TarArchiveInputStream(in)
    private var fileReader: Option[BufferedScanner] = None
    private var _lineCount = 0

    private def getReader: Option[BufferedScanner] = {
      def nextReader: Option[BufferedScanner] = {
        val entry = tarIn.getNextTarEntry
        if (entry == null)
          None
        else if (entry.isDirectory)
          nextReader
        else if (isFASTAFile(entry.getName)) {
          fileReader = Some(new BufferedScanner(tarIn, bufferSize))
          fileReader
        }
        else
          nextReader
      }

      fileReader.flatMap {
        reader =>
          if (reader.LA(1) == BufferedScanner.EOF) {
            // End of the file
            fileReader = None
          }
          fileReader
      }.orElse(nextReader)
    }

    def nextLine: CharSequence = {
      getReader.map {
        reader =>
          _lineCount += 1
          reader.nextLine()
      }.getOrElse(null)
    }

    def LA1 = getReader.map(_.LA(1)).getOrElse(BufferedScanner.EOF)

    def lineCount = _lineCount

    def close {
      tarIn.close
      in.close
    }
  }


  private def createStream(reader: LineReader): Stream[FASTASequenceReader] = {

    def loop(prevStream: Option[SequenceStream]): Stream[FASTASequenceReader] = {
      prevStream.foreach(_.invalidate)
      def createStream: Stream[FASTASequenceReader] = {
        val next = reader.LA1
        next match {
          case BufferedScanner.EOF => Stream.empty[FASTASequenceReader]
          case '>' => {
            // Start of the description line
            val ss = new SequenceStream(reader)
            val fr = new FASTASequenceReader(reader.nextLine.toString, ss)
            fr #:: loop(Some(ss))
          }
          case '#' => {
            // Skip the comment line
            reader.nextLine
            createStream
          }
          case _ => {
            val line = reader.nextLine
            warn("invalid data at line %,d:\n", reader.lineCount, line)
            createStream
          }
        }
      }
      createStream
    }
    loop(None)
  }


  /**
   * Open and read the fasta sequences from file.
   * The stream passed to the function f is not thread-safe because the stream is dynamically created by reading the file from head to tail.
   *
   * @param fileName
   * @return
   */
  def read[A](fileName: String)(f: Stream[FASTASequenceReader] => A): A = {
    read(createLineReader(fileName))(f)
  }

  def read[A](fastaData: Reader)(f: Stream[FASTASequenceReader] => A): A = {
    val r = new StandardLineReader(new BufferedScanner(fastaData, 4 * 1024 * 1024))
    read(r)(f)
  }

  def read[A](fastaData: InputStream)(f: Stream[FASTASequenceReader] => A): A = {
    val r = new StandardLineReader(new BufferedScanner(fastaData, 4 * 1024 * 1024))
    read(r)(f)
  }

  private def read[A](lineReader: LineReader)(f: Stream[FASTASequenceReader] => A): A = {
    val stream = createStream(lineReader)
    try {
      f(stream)
    }
    finally {
      lineReader.close()
    }
  }

  def isFASTAFile(fileName: String) = hasExt(fileName, ".fa", ".fasta", ".fs")

  private def hasExt(fileName: String, extList: String*) = extList.exists(fileName.endsWith(_))

  private def createLineReader(fileName: String, bufferSize:Int = DEFAULT_BUFFER_SIZE): LineReader = {
    if (hasExt(fileName, ".tgz", ".tar.gz"))
      new TarFileLineReader(new BufferedInputStream(new GZIPInputStream(new FileInputStream(fileName)), bufferSize), bufferSize)
    else if (hasExt(fileName, ".tar"))
      new TarFileLineReader(new BufferedInputStream(new FileInputStream(fileName), bufferSize), bufferSize)
    else
      new StandardLineReader(new BufferedScanner(new FileInputStream(fileName), bufferSize))
  }


}


/**
 * Reader of each FASTA entry
 * @param description
 * @param ss
 */
class FASTASequenceReader(val description: String, private val ss: FASTA.SequenceStream) {
  import FASTA._
  /**
   * name of the sequence
   */
  lazy val name = extractSequenceNameFrom(description)

  /**
   * Stream for reading genome sequences line by line. This stream can be used only once.
   */
  val stream: Stream[String] = ss.toStream

  /**
   * Create the String representation of the sequence
   */
  lazy val sequence: String = {
    if (!ss.isValidStream)
      sys.error("This sequence was already read somewhere else")

    val b = new StringBuilder
    for (line <- stream) {
      b.append(line.trim)
    }
    b.result
  }
}
