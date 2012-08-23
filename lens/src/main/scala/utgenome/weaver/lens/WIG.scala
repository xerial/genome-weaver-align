package utgenome.weaver.lens

import util.parsing.combinator.RegexParsers
import java.io.{FileReader, BufferedReader}
import xerial.core.log.Logging

//--------------------------------------
//
// WIG.scala
// Since: 2012/07/18 12:49
//
//--------------------------------------

/**
 * WIG format
 * @author leo
 */
object WIG {

  sealed abstract class Header extends WIG

  case class Browser(line: String) extends Header
  case class Track(property: Map[String, String]) extends Header
  case class Comment(line: String) extends WIG
  case class VariableStep(chrom: String, span: Int = 1) extends Header
  case class FixedStep(chrom: String, start: Int, step: Int = 1, span: Int = 1) extends Header

  sealed abstract class Data extends WIG
  case class VariableStepValue(position: Int, value: Float) extends Data
  case class FixedStepValue(value: Float) extends Data

  case class Error(message:String) extends WIG
}

sealed abstract class WIG

/**
 * WIG format parser.
 *
 * Wiggle format (WIG)
 * http://genome.ucsc.edu/goldenPath/help/wiggle.html
 */
object WIGParser extends RegexParsers with Logging {

  /**
   * remove quotation symbols
   * @param s
   * @return
   */
  protected def unquote(s: String): String = s.substring(1, s.length() - 1)

  def paramName: Parser[String] = """[A-Za-z0-9:\-_\.]+""".r

  def stringLiteral: Parser[String] =
    ("\"" + """([^"\p{Cntrl}\\]|\\[\\/bfnrt]|\\u[a-fA-F0-9]{4})*""" + "\"").r ^^ {
      unquote(_)
    }

  def quotation: Parser[String] =
    ("'" + """([^'\p{Cntrl}\\]|\\[\\/bfnrt]|\\u[a-fA-F0-9]{4})*""" + "'").r ^^ {
      unquote(_)
    }

  def value: Parser[String] = """[^\"'\s]+""".r

  def paramValue: Parser[String] = stringLiteral | quotation | value

  def param: Parser[(String, String)] = paramName ~ "=" ~ paramValue ^^ {
    case key ~ "=" ~ value => (key, value)
  }

  def header: Parser[(String, Map[String, String])] = paramName ~ rep(param) ^^ {
    case p ~ params => (p, Map() ++ params)
  }

  def parseHeader(line: String) : Either[NoSuccess, (String, Map[String, String])] = {
    parseAll(header, line) match {
      case Success(result, next) => Right(result)
      case failure : NoSuccess => Left(failure)
    }
  }

  class LineIterator(in:BufferedReader) extends Iterator[String] {
    private var nextLine : String = null
    def hasNext = {
      if(nextLine == null)
        nextLine = in.readLine
      nextLine != null
    }
    def next : String = {
      if(hasNext) {
        val line = nextLine
        nextLine = null
        line
      }
      else
        Iterator.empty.next
    }
  }


  def open[U](fileName:String)(body:Iterator[String] => U) {
    val in = new BufferedReader(new FileReader(fileName))
    try
      body(new LineIterator(in))
    finally
      in.close
  }


  def parse(file:String) {
    open(file) { f =>
      for(line <- f) {
        val wig = parseLine(line)
        trace(wig)
      }
    }
  }

  def parseLine(line: String): WIG = {

    def convert[A](s:String, f:String => A) : Either[Throwable, A] = 
      scala.util.control.Exception.allCatch.either(f(s))
    
    def toInt(s:String) = convert(s, {_.toInt})
    def toFloat(s:String) = convert(s, {_.toFloat})

    def parse : Either[NoSuccess, WIG] = {
      if (line.startsWith("#"))
        Right(WIG.Comment(line))
      else if (line.startsWith("browser"))
        Right(WIG.Browser(line))
      else if (line.startsWith("track"))
        parseHeader(line).right.map(header => WIG.Track(header._2))
      else if (line.startsWith("variableStep")) {
        parseHeader(line).right.map{
          case (name, props) =>
            (props.get("chrom"), toInt(props.getOrElse("span", "1"))) match {
              case (Some(chr), Right(sp)) => WIG.VariableStep(chrom=chr, span=sp)
              case _ => WIG.Error("invalid line: " + line)
            }
        }
      }
      else if (line.startsWith("fixedStep")) {
        parseHeader(line).right.map{
          case (name, props) =>
            val chrom = props.get("chrom")
            val start = props.get("start")
            val step = toInt(props.get("step").getOrElse("1"))
            val span = toInt(props.get("span").getOrElse("1"))
            (chrom, start, step, span) match {
              case (Some(chr), Some(s), Right(st), Right(sp)) =>
                WIG.FixedStep(chrom=chr, start=s.toInt, step=st, span=sp)
              case _ => WIG.Error("invalid line: " + line)
            }
        }
      }
      else {
        // data line
        def invalidLine = WIG.Error("invalid line: " + line)
        val c = line.trim.split("""\s+""")
        val r = c match {
          case Array(step, value) => (toInt(step), toFloat(value)) match {
            case (Right(st), Right(v)) => WIG.VariableStepValue(st, v)
            case _ => invalidLine
          }
          case Array(value) => toFloat(value) match {
            case Right(v) => WIG.FixedStepValue(value.toFloat)
            case _ => invalidLine
          }
          case _ => invalidLine
        }
        Right(r)
      }
    }

    parse match {
      case Right(m) => m
      case Left(error) => WIG.Error(error.toString)
    }
  }

}


