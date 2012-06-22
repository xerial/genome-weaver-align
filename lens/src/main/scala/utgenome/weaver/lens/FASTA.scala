//--------------------------------------
//
// FASTA.scala
// Since: 2012/06/22 10:00 AM
//
//--------------------------------------

package utgenome.weaver.lens

import xerial.silk.util.io.{BufferedScanner, FileSource}
import java.io._
import scala.Some
import annotation.tailrec
import xerial.silk.util.Logger
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import java.util.zip.GZIPInputStream


/**
 * FASTA file reader.
 *
 * <h3>Usage</h3>
 * <pre>
 * <code>
 * FASTA.read(randomFASTAReader) { stream =>
 *   // reader for each chromosome
 *   for(r : FASTAEntryReader <- stream) {
 *      val name = r.name
 *      val description = r.description
 *      // read each sequence line
 *      for(line <- r.lines)
 *          ...
 *   }
 * </code>
 * </pre>
 *
 * If you need to retrieve the entire sequence as a single String, use [[utgenome.weaver.lens.FASTAEntryReader#sequence]].
 *
 * @author leo
 */
object FASTA extends Logger {


  private def create2bitIndexFrom(stream:Stream[FASTAEntryReader]) : FASTAIndex2bit = {
    val index = Array.newBuilder[FASTAEntryIndex]
    var offset = 0L
    val b = ACGTSeq.newBuilder
    for(r : FASTAEntryReader <- stream) {
      val desc = r.description
      val name = r.name
      debug("loading %s", name)
      for(line <- r.lines) {
        b += line
      }
      // TODO preserve bit vector of Ns
      index += new FASTAEntryIndex(name, desc, offset, (b.numBases - offset).toInt)
      offset = b.numBases
    }
    new FASTAIndex2bit(b.result, index.result)
  }

  def create2bitIndexFrom(fastaFile:String) : FASTAIndex2bit = read(fastaFile)(create2bitIndexFrom)
  def create2bitIndexFromTarGZ(in:InputStream) = readTarGZ(in)(create2bitIndexFrom)
  def create2bitIndexFromTarGZ(in:Reader) = read(in)(create2bitIndexFrom)


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

  // 4MB
  val DEFAULT_BUFFER_SIZE = 4 * 1024 * 1024



  private[glens] class TarFileLineReader(in: InputStream, bufferSize: Int = DEFAULT_BUFFER_SIZE) extends LineReader {
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


  private def createStream(reader: LineReader): Stream[FASTAEntryReader] = {

    def loop(prevStream: Option[SequenceStream]): Stream[FASTAEntryReader] = {
      prevStream.foreach(_.invalidate)
      def createStream: Stream[FASTAEntryReader] = {
        val next = reader.LA1
        next match {
          case BufferedScanner.EOF => Stream.empty[FASTAEntryReader]
          case '>' => {
            // Start of the description line
            val ss = new SequenceStream(reader)
            val fr = new FASTAEntryReader(reader.nextLine.toString, ss)
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
   * Open and read the fasta sequences from file. This method supports plain fasta (.fasta, .fa, .fs) and compressed (.tar.gz, .gz) files.
   *
   * The stream passed to the function f is not thread-safe because the stream is dynamically created by reading the file from head to tail.
   *
   * @param fileName
   * @return
   */
  def read[A](fileName: String)(f: Stream[FASTAEntryReader] => A): A = read(fileName, DEFAULT_BUFFER_SIZE)(f)
  def read[A](fileName: String, bufferSize: Int)(f: Stream[FASTAEntryReader] => A): A = {
    read(createLineReader(fileName, bufferSize))(f)
  }

  def read[A](fastaData: Reader)(f: Stream[FASTAEntryReader] => A): A = read(fastaData, DEFAULT_BUFFER_SIZE)(f)
  def read[A](fastaData: Reader, bufferSize: Int)(f: Stream[FASTAEntryReader] => A): A = {
    val r = new StandardLineReader(new BufferedScanner(fastaData, bufferSize))
    read(r)(f)
  }


  def read[A](fastaData: InputStream)(f: Stream[FASTAEntryReader] => A): A = read(fastaData, DEFAULT_BUFFER_SIZE)(f)
  def read[A](fastaData: InputStream, bufferSize: Int)(f: Stream[FASTAEntryReader] => A): A = {
    val r = new StandardLineReader(new BufferedScanner(fastaData, 4 * 1024 * 1024))
    read(r)(f)
  }



  def readTarGZ[A](fastaTGZData: InputStream)(f: Stream[FASTAEntryReader] => A): A = readTarGZ(fastaTGZData, DEFAULT_BUFFER_SIZE)(f)
  def readTarGZ[A](fastaTGZData: InputStream, bufferSize: Int)(f: Stream[FASTAEntryReader] => A): A = {
    val r = new TarFileLineReader(new BufferedInputStream(new GZIPInputStream(fastaTGZData), bufferSize), bufferSize)
    read(r)(f)
  }

  private[glens] def read[A](lineReader: LineReader)(f: Stream[FASTAEntryReader] => A): A = {
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

  private def createLineReader(fileName: String, bufferSize: Int = DEFAULT_BUFFER_SIZE): LineReader = {
    if (hasExt(fileName, ".tgz", ".tar.gz"))
      new TarFileLineReader(new BufferedInputStream(new GZIPInputStream(new FileInputStream(fileName)), bufferSize), bufferSize)
    else if (hasExt(fileName, ".tar"))
      new TarFileLineReader(new BufferedInputStream(new FileInputStream(fileName), bufferSize), bufferSize)
    else
      new StandardLineReader(new BufferedScanner(new FileInputStream(fileName), bufferSize))
  }


}


/**
 * Reader of each FASTA entry. Note that [[utgenome.weaver.lens.FASTAEntryReader#lines]] is evaluated lazily.
 * @param description
 * @param ss
 */
class FASTAEntryReader(val description: String, private val ss: FASTA.SequenceStream) {

  import FASTA._

  /**
   * name of the sequence
   */
  val name = extractSequenceNameFrom(description)

  /**
   * Get a stream for reading the genome sequence line by line.
   * This stream can be used only once.
   */
  val lines: Stream[String] = ss.toStream

  /**
   * Extract the entire sequence. This result can be large. If you want to process genome sequences line by line, use [[utgenome.weaver.lens.FASTAEntryReader#lines]] method.
   */
  lazy val sequence: String = {
    if (!ss.isValidStream)
      sys.error("This sequence was already read somewhere else")

    val b = new StringBuilder
    for (line <- lines) {
      b.append(line.trim)
    }
    b.result
  }
}


