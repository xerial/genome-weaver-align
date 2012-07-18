package xerial.silk.glens

import xerial.silk.util.SilkSpec

//--------------------------------------
//
// WIGTest.scala
// Since: 2012/07/18 13:49
//
//--------------------------------------

/**
 * @author leo
 */
class WIGTest extends SilkSpec {

  import WIGParser._

  def parse(line:String) : WIG = {
    val r = parseLine(line)
    debug(r)
    r
  }

  "WIG" should {

    "parse param" in {
      val r = WIGParser.parseAll(paramName, "chrom")
      debug(r)
    }

    "parse value" in {
      val r = WIGParser.parseAll(paramValue, "chr19")
      debug(r)

      val r2 = parseAll(paramValue, """'hello world'""")
      debug(r2)
    }
    
    "parse key value pair" in {
      debug("parse key-value")
      val r = WIGParser.parseAll(param, "chrom=chr19")
      debug(r)
    }
    
    "parse variableStep line" in {
      debug("parser variableStep")
      val r = WIGParser.parseAll(header, "variableStep chrom=chr19 span=150")
      debug(r)
    }

    "parse fixedStep line" in {

      val r = WIGParser.parseHeader("fixedStep chrom=19 start=10000 step=300 span=200")
      debug(r)
    }

    "parse lines" in {
      parse("track type=wiggle_0 name=\"fixedStep\" description=\"fixedStep format\" visibility=full autoScale=off viewLimits=0:1000 color=0,200,100 maxHeightPixels=100:50:20 graphType=points priority=20")
      parse("browser position chr19:49304200-49310700")
      parse("browser hide all")

      parse("variableStep chrom=chr19 span=150")
      parse("fixedStep chrom=chr19 start=49307401 step=300 span=200")
      parse("49304701 10.0")
      parse("1000")
    }

    "parse lines with errors" in {
      parse("track type=wiggle_0 name=fa09r3 43")
    }

  }

}