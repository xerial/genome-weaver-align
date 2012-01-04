/*
 * Copyright 2012 Taro L. Saito
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.utgenome.weaver
import org.xerial.util.log.{Logger => XerialLogger}
import collection.mutable.WeakHashMap

//--------------------------------------
//
// Logger.scala
// Since: 2012/01/04 13:17
//
//--------------------------------------

object LoggerHolder {
  private val holder = new WeakHashMap[Class[_], XerialLogger]
  def logger(x:Any) : XerialLogger = {
    val c = x.getClass
    if(holder.contains(c))
      holder(c)
    else {
      val l = XerialLogger.getLogger(c)
      holder += (c -> l)
      l
    }
  }
}

/**
 * Enabling logging in the class that integrates this trait
 *
 *
 * Created at 2012/01/04 13:17
 * @author leo
 */
trait Logger {
  private[this] lazy val _logger = LoggerHolder.logger(this)

  private def wrap(format: => String, args: Any*)(tester: => Boolean, logMethod: String => Boolean) : Boolean = {
     if(tester) {
       logMethod {
         if (args.size == 0)
           format
         else  {
           format.format(args:_*)
         }
       }
     }
     else false
  }
  
  type LogFunction = (=>String, Any*) => Boolean

  def fatal(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isFatalEnalbed, _logger.fatal)
  def error(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isErrorEnabled, _logger.error)
  def warn(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isWarnEnabled, _logger.warn)
  def info(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isInfoEnabled, _logger.info)
  def debug(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isDebugEnabled, _logger.debug)
  def trace(format: => String, args: Any*) = wrap(format, args:_*)(_logger.isTraceEnabled, _logger.trace)


}