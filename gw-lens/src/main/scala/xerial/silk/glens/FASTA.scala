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


/**
 * FASTA file reader
 *
 * @author leo
 */
object FASTA {

  def extractSequenceNameFrom(descriptionLine: String) = {
    val trimmed = descriptionLine.substring(1).dropWhile(c => Character.isWhitespace(c))
    trimmed.takeWhile(c => !Character.isWhitespace(c) && c != '|')
  }

  /**
   * Reading
   * @param in
   */
  private[glens] class SequenceStream(in: BufferedScanner) {
    private var finishedReading = false
    private var initialized = false

    private def noMoreData: Boolean = {
      val ch = in.LA(1)
      if (!finishedReading && (ch == '>' || ch == BufferedScanner.EOF))
        finishedReading = true
      finishedReading
    }

    def isValidStream = !initialized && !finishedReading

    private[FASTA] def invalidate {
      while (!noMoreData)
        in.nextLine()
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
          in.nextLine().toString #:: loop
      loop
    }
  }


  private def createStream(in: BufferedScanner): Stream[FASTASequenceReader] = {
    def loop(prevStream: Option[SequenceStream]): Stream[FASTASequenceReader] = {
      prevStream.foreach(_.invalidate)
      def createStream: Stream[FASTASequenceReader] = {
        val next = in.LA(1)
        next match {
          case BufferedScanner.EOF => Stream.empty[FASTASequenceReader]
          case '>' => {
            // Start of the description line
            val ss = new SequenceStream(in)
            val s = new FASTASequenceReader(in.nextLine.toString, ss)
            s #:: loop(Some(ss))
          }
          case '#' => {
            // Skip the comment line
            in.nextLine
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
  def read[A](fileName: String)(f: Stream[FASTASequenceReader] => A): A =
    read(new FileInputStream(fileName))(f)


  def read[A](fastaData: InputStream)(f: Stream[FASTASequenceReader] => A): A = {
    val scanner = new BufferedScanner(fastaData, 4 * 1024 * 1024)
    val stream = createStream(scanner)
    try {
      f(stream)
    }
    finally {
      scanner.close()
    }
  }


}


/**
 * Reader of each FASTA entry
 * @param description
 * @param ss
 */
class FASTASequenceReader(val description: String, private val ss: FASTA.SequenceStream) {
  /**
   * name of the sequence
   */
  lazy val name = FASTA.extractSequenceNameFrom(description)

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
