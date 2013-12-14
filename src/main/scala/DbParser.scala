package org.hazy.gibbs

import scala.collection.immutable.{ Map, List }
import scala.collection.mutable._
import scala.util.Random
import scala.math.{ abs, exp, log }
import scala.io.Source
import java.io._

object DBParser {

  // *****************************************
  // * Helper functions to Parse DB's format *
  // *****************************************
  def parseWeight(szWeightName: String): Map[Int, Weight] = {
    return Source.fromFile(szWeightName).getLines().map(_.split("\t")).map(
      v =>
        // [WEIGHT ID] [VALUE] [IS_FIXED_WEIGHT]
        // new Weight(v(0).toInt, v(1).toDouble, v(2).toBoolean)
        (v(0).toInt, new Weight(v(0).toInt, v(1).toDouble))).toMap
  }

  def parseVariables(szVarName: String, szFactorName: String,
    htW: Map[Int, Weight]): FactorManager = {
    val ht = new HashMap[Int, Set[(Int, Variable, Int)]] 
    val vs = new HashMap[Int, Variable] 
    val fm = new FactorManager()
    // 0[VARIABLE_ID] 1[FACTOR_ID] 2[POSITION] 3[IS_POSITIVE] 4[DATA_TYPE] 5[INIT_VALUE] 6[IS_EVIDENCE] 7[IS_QUERY]
    val iter = Source.fromFile(szVarName).getLines().map(_.split("\t"))
    for (v <- iter) {
      var varId: Int = v(0).toInt
      val value = (v(5) != 0)
      var x = new Variable(varId, value)
      vs(varId) = x

        val fid = v(1).toInt
        val pos = v(2).toInt
        val sign = if (v(3).toBoolean) 1 else 0
        Util.checkAndAddSet[Int, (Int, Variable, Int)](ht, fid, (pos, x, sign))
    }
    println("Variables Created.")

    // Now we can parse the Factors
    // [FACTOR_ID] [WEIGHT ID] [FUNCTION_TYPE]
    val factIt = Source.fromFile(szFactorName).getLines().map(_.split("\t"))
    for (f <- factIt) {
      val fid = f(0).toInt
      val funId = f(2)
      val wid = f(1).toInt

      val args = ht(fid).toList.sortWith(_._1 < _._1)
      //val signMap = args.map( q => (q._2, (q._3 < 0))).toMap
      fm.addFactor(new ImplyFactor(fid, args.map(_._2).toArray, htW(wid)))
    }
    return fm
  }

  // Main parsing function for DB's format
  def parse(szWeightName: String, szVarName: String, szFactorName: String): FactorManager = {
    val htW = parseWeight(szWeightName)
    return parseVariables(szVarName, szFactorName, htW)
  }
}