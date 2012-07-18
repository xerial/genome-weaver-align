package xerial.silk.glens

import util.parsing.combinator.RegexParsers
import xerial.silk.util.Logger

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
object WIGParser extends RegexParsers with Logger {

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

  def parseHeader(line: String) : Either[(String, Map[String, String]), NoSuccess] = {
    parseAll(header, line) match {
      case Success(result, next) => Left(result)
      case failure : NoSuccess => Right(failure)
    }
  }


  def parseLine(line: String): WIG = {

    def parse : Either[WIG, NoSuccess] = {
      if (line.startsWith("#"))
        Left(WIG.Comment(line)) // comment
      else if (line.startsWith("browser"))
        Left(WIG.Browser(line))
      else if (line.startsWith("track"))
        parseHeader(line).left.map(header => WIG.Track(header._2))
      else if (line.startsWith("variableStep")) {
        parseHeader(line).left.map{
          case (name, props) =>
            val v = props.get("chrom").map{chrom =>
              val span = props.getOrElse("span", "1").toInt
              WIG.VariableStep(chrom=chrom, span=span)
            }
            v.getOrElse(WIG.Error("missing column in: " + line))
        }
      }
      else if (line.startsWith("fixedStep")) {
        parseHeader(line).left.map{
          case (name, props) =>
            val chrom = props.get("chrom")
            val start = props.get("start")
            val step = props.get("step").getOrElse("1").toInt
            val span = props.get("span").getOrElse("1").toInt

            (chrom, start, step, span) match {
              case (Some(chr), Some(s), st:Int, sp:Int) =>
                WIG.FixedStep(chrom=chr, start=s.toInt, step=st, span=sp)
              case _ => WIG.Error("missing value in: " + line)
            }
        }
      }
      else {
        // data line
        val c = line.split("""\s+""")
        c match {
          case Array(step, value) => Left(WIG.VariableStepValue(step.toInt, value.toFloat))
          case Array(value) => Left(WIG.FixedStepValue(value.toFloat))
          case _ => Left(WIG.Error("invalid line: " + line))
        }
      }
    }

    parse match {
      case Left(m) => m
      case Right(error) => WIG.Error(error.toString)
    }
  }

}


