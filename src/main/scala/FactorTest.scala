package org.hazy.gibbs

import scala.collection.immutable.{ Map, List }
import scala.collection.mutable._
import scala.util.Random
import scala.math.{ abs, exp, log }
import scala.io.Source
import java.io._

//
// Utils
//
object Util {
  def checkAndAddSet[K, V](ht: HashMap[K, Set[V]], k: K, v: V) {
    if (!(ht contains k)) { ht(k) = new HashSet[V]() }
    ht(k).add(v)
  }

  def printToFile(f: java.io.File)(op: java.io.PrintWriter => Unit) {
    val p = new java.io.PrintWriter(f)
    try { op(p) } finally { p.close() }
  }
}
// ***************************************************
// There are three major parts of factor graphs
// Variables, Factors, and Weights.
// These are the classes.

// TODO: Template the variable types
class Variable(id0: Int, value0: Boolean) {
  var _value = value0
  def value: Boolean = return _value
  def id = id0
  //def setTrue() = synchronized { _value = true }
  //def setFalse() = synchronized { _value = false }
  def setTrue() = { _value = true }
  def setFalse() = { _value = false }

  override def toString() =
    "" + id + "->" + _value
}

//
// Weights may be tied together across several factors.
//
class Weight(id: Int, w0: Double) {
  def id0 = id
  var w = w0 // this is a log-scale weight
}

// 
// the factor returns the weight if any variable is true 
// 
abstract class Factor(id0: Int, vars0: Array[Variable], w0: Weight) {
  def id = id0
  def vars = vars0
  def w = w0 // this could be tied
  /* Note variable weights are different from weights */
  def WeightVariable(u: Variable): (Double, Double)
  def EvaluateFactor(): Boolean /* evaluate this factor */
}

// 
// Several simple types of Factors
//
class OrFactor(id0 : Int, vars0: Array[Variable], w0 : Weight) extends Factor(id0 : Int, vars0: Array[Variable], w0 : Weight) { 
  override def WeightVariable(u : Variable) : (Double,Double) = { 
    // If everythign aside from u is false, then u is influential
    if (vars.filter(x => x != u).forall(v => !v.value)) { 
      return (w.w,0.0)
    }
    return (0.0,0.0)
  }
  override def EvaluateFactor() : Boolean = { return vars.exists(v => v.value) }
}

class AndFactor(id0 : Int, vars0: Array[Variable], w0 : Weight) extends Factor(id0 : Int, vars0: Array[Variable], w0 : Weight) { 
  override def WeightVariable(u : Variable) : (Double,Double) = { 
    // if everything but u is true, then u is influential.
    if (vars.filter(x => x != u).forall(v => v.value)) { 
      return (w.w,0.0)
    }
    return (0.0,0.0)
  }
  override def EvaluateFactor() : Boolean = { return vars.forall(v => v.value) }
}

class ImplyFactor(id0 : Int, vars0: Array[Variable], w0 : Weight) extends Factor(id0 : Int, vars0: Array[Variable], w0 : Weight) {
  var imps   = vars0.slice(0,vars.size - 1)
  var output = vars0(vars.size - 1)

  override def WeightVariable(u : Variable) : (Double,Double) = { 
    // if either the body is false or the output ...
    if (!imps.filter(x => x != u).forall(v => v.value) || ((output != u) && output.value)) { 
      return (0.0,0.0)
    } else { 
      if( output == u ) { 
	return (w.w,0.0)
      } else { 
	return (0.0,w.w)
      }
    }
  } 
  override def EvaluateFactor() : Boolean = { 
    val n = vars.size
    return vars.take(n - 1).forall(v => !v.value) || vars(n-1).value
  }
}


// ****************************************************
// A Ce factor subsumes the above hardcoded factors. 
// A Ce factor is an Orfactor but allows Not gates (signs) for each input.
class CeFactor(id0: Int, vars0: Array[Variable], w0: Weight, signs0: Map[Variable, Boolean]) extends Factor(id0: Int, vars0: Array[Variable], w0: Weight) {
  def signs = signs0

  override def WeightVariable(u: Variable): (Double, Double) = {
    //synchronized {
      // If everything aside from u is false (after signing), then u is influential
      if (vars.filter(x => x != u).forall(v => !(v.value ^ signs(v)))) {
        return if (signs(u)) (0.0, w.w) else (w.w, 0.0)
      }
      return (0.0, 0.0)
    //}
  }
  override def EvaluateFactor(): Boolean = { return vars.exists(v => v.value ^ signs(v)) }
}

// 
// A FactorManager keeps track of the associations between variables, factors, and weights.
// It is also where we put in samplers for now
class FactorManager {
  val factorMap = new HashMap[Variable, Set[Factor]]
  val weightMap = new HashMap[Weight, Set[Factor]]
  val variables = new HashSet[Variable]
  def r = Random

  def registerVar(v: Variable, f: Factor) {
    if (!(variables contains v)) {
      variables.add(v)
      factorMap(v) = new HashSet[Factor]()
    }
    factorMap(v).add(f)
  }

  def registerWeight(f: Factor) {
    if (!(weightMap contains f.w)) {
      weightMap(f.w) = new HashSet[Factor]()
    }
    weightMap(f.w).add(f)
  }

  def addFactor(f: Factor) {
    for (v <- f.vars) {
      registerVar(v, f)
      registerWeight(f)
    }
  }

  // This is the individual uncollapsed Gibbs sampler code.
  def sampleVariable(v: Variable) {
    var (pos: Set[Double], neg: Set[Double]) = factorMap(v).map(_.WeightVariable(v)).unzip
    var pos0 = pos.foldLeft(0.0)(_ + _)
    var neg0 = neg.foldLeft(0.0)(_ + _)
    var s = r.nextDouble()

    // if( s*(exp(pos0)+exp(neg0)) <= exp(pos0) ) { v.setTrue() } else { v.setFalse() }
    if (s * (1 + exp(neg0 - pos0)) <= 1.0) { v.setTrue() } else { v.setFalse() }
  }

  // Sample all known variables
  def sampleVariables() { for (v <- r.shuffle(variables.toList)) { sampleVariable(v) } }

  // Sample variables that are in the mask
  def sampleVariablesMask(mask: Set[Variable]) {
    for (v <- r.shuffle(variables.toList.filter(u => mask contains u))) {
      sampleVariable(v)
    }
  }

  // Parallel sampling
  def sampleVariablesMaskParallel(masks: Array[MaskVariable], nThreads: Int) {
    var t = new Array[Thread](nThreads)
    // For all threads, set them up
    for (i <- 0 until nThreads) {
      t(i) = new Thread(new Sample(this, masks(i).variables, i))
      t(i).start()
    }

    // Joins all the threads
    for (i <- 0 until nThreads) { t(i).join() }
    
    // Merges the result
    for (i <- 0 until nThreads) {
      t(i) = new Thread(new MergeResult(this, masks(i).variables, i))
      t(i).start()
    }
    
    // Joins all the threads
    for (i <- 0 until nThreads) { t(i).join() }

  }

  /* Calculates the probability that a variable is true */
  def marginals(nSamples: Int, mask: Set[Variable]): List[(Variable, Double)] = {
    var trueCount = new HashMap[Variable, Int]()
    for (v <- variables) { trueCount(v) = 0 }

    for (i <- 0 until nSamples) {
      sampleVariablesMask(mask)

      for (v <- variables) {
        if (v.value) { trueCount(v) += 1 }
      }
    }
    return trueCount.toList.map(x => (x._1, x._2 / nSamples.toDouble))
  }

  //
  // this computes the marginal probability that each factor is true
  // 
  def marginalsFactor(nSamples: Int, mask: Set[Variable], factorSet: scala.collection.immutable.Set[Factor]): Map[Factor, Double] = {
    var trueCount = new HashMap[Factor, Int]()
    for (f <- factorSet) { trueCount(f) = 0 }

    for (i <- 0 until nSamples) {
      sampleVariablesMask(mask)

      for (f <- factorSet) {
        if (f.EvaluateFactor()) { trueCount(f) += 1 }
      }
    }
    return trueCount.toList.map(x => (x._1, x._2 / nSamples.toDouble)).toMap
  }

  def histogram(nSamples: Int, mask: Set[Variable]) {
    var trueCount = new HashMap[Variable, Int]()
    for (v <- variables) { trueCount(v) = 0 }

    for (i <- 0 until nSamples) {
      sampleVariablesMask(mask)

      for (v <- variables) {
        if (v.value) { trueCount(v) += 1 }
      }
    }
    println(trueCount.toList.map(x => (x._1.id, x._2)) + " of " + nSamples)
  }

  // histogram file
  def MarginalsFile(nSamples: Int, mask: Set[Variable], nThreads: Int, szOut: String) {
    var trueCount = new HashMap[Variable, Int]()
    val startTime = System.nanoTime
    for (v <- variables) { trueCount(v) = 0 }

    //val nThreads = 1

    var masks = new Array[MaskVariable](nThreads)
    for (i <- 0 until nThreads) masks(i) = new MaskVariable(Set.empty[Variable])

    // Partition variables into different masks
    var i = 0
    for (v <- r.shuffle(variables.toList.filter(mask contains _))) {
      // allocate v to i % nThreads mask
      masks(i % nThreads).variables.add(v)
      i += 1
    }
    /*val v1 = r.shuffle(variables.toList.filter(mask contains _))
    var it = v1.grouped(v1.length / nThreads);
    for (i <- 0 until nThreads) {
      masks(i).variables ++= it.next
    }*/
    //for (i <- 0 until nThreads) masks(i).variables ++= mask
    val time1 = System.nanoTime
    println("partition " + (time1 - startTime) / 1e9)

    for (i <- 0 until nSamples) {
      //var time2 = System.nanoTime
      sampleVariablesMaskParallel(masks, nThreads)

      //var time3 = System.nanoTime
      //println("sample " + (time3 - time2) / 1e9)
      
      //mask.clear()
      //for (i <- 0 until nThreads) mask ++= masks(i).variables
      
      for (v <- variables) {
        if (v.value) { trueCount(v) += 1 }
      }
      if (i % 10 == 1) { print(".") }
    }
    val endTime = System.nanoTime
  
    // println(trueCount.toList.map( x => (x._1.id,x._2)) + " of " + nSamples)
    // val data = trueCount.toList.map( x => (x._1.id,x._2, x._2/nSamples.toDouble))
    val time_s = (endTime - startTime) / 1e9
    println("\ntime = " + time_s + " iterations = " + nSamples + " active vars = " + mask.size)
    println("samples/s = " + nSamples / time_s + " tokens/s = " + (nSamples * mask.size) / time_s)

    val cesFormat = trueCount.toList.map(x => x._1.id + "\t" + 1 + "\t" + (x._2 / nSamples.toDouble))
    Util.printToFile(new File(szOut))(p => {
      cesFormat.foreach(p.println)
    })
  }

}

class Sample(fm: FactorManager, mask: Set[Variable], id0: Int) extends Runnable {
  def id = id0
  def r = Random

  def run() {
    for (v <- mask) {
      var (pos: Set[Double], neg: Set[Double]) = fm.factorMap(v).map(_.WeightVariable(v)).unzip
      var pos0 = pos.foldLeft(0.0)(_ + _)
      var neg0 = neg.foldLeft(0.0)(_ + _)
      var s = r.nextDouble()
      if (s * (1 + exp(neg0 - pos0)) <= 1.0) { v.setTrue() } else { v.setFalse() }
    }
    Thread.`yield`()
  }
}

class MergeResult(fm: FactorManager, mask: Set[Variable], id0:Int) extends Runnable {
  def run() {
    fm.variables ++= mask
    Thread.`yield`()
  }
}

class MaskVariable(v0: Set[Variable]) {
  def variables = v0
}

object CeParser {

  // *****************************************
  // * Helper functions to Parse Ce's format *
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
    val ht = new HashMap[Int, Set[(Int, Variable, Int)]] /* what's this? */
    val vs = new HashMap[Int, Variable] /* variables */
    val fm = new FactorManager()
    // [VARIABLE_ID] [PROPERTY] [LOWER_BOUND] [UPPER_BOUND] [NUMBER_OF_FACTOR] [FACTOR_INFO …] [INIT VALUE]
    val iter = Source.fromFile(szVarName).getLines().map(_.split("\t"))
    for (v <- iter) {
      var varId: Int = v(0).toInt
      val prop = v(1)
      val nFactors = v(4).toInt
      val value = (v(4 * nFactors + 5) != 0)
      var x = new Variable(varId, value)
      vs(varId) = x

      for (j <- 0 until nFactors) {
        //[FACTOR ID] [GROUP_ID] [POSITION] [SIGN]
        val fid = v(4 * j + 5).toInt
        val pos = v(4 * j + 5 + 2).toInt
        val sign = v(4 * j + 5 + 3).toInt
        Util.checkAndAddSet[Int, (Int, Variable, Int)](ht, fid, (pos, x, sign))
      }
    }
    println("Variables Created.")

    // Now we can parse the Factors
    // [FACTOR_ID] [NUM_VARIABLE_IN_FACTOR] [FACTOR_FUN_ID] [WEIGHT ID]
    val factIt = Source.fromFile(szFactorName).getLines().map(_.split("\t"))
    for (f <- factIt) {
      val fid = f(0).toInt
      val nFacts = f(1).toInt
      val funId = f(2)
      val wid = f(3).toInt

      val args = ht(fid).toList.sortWith(_._1 < _._1)
      val signMap = args.map( q => (q._2, (q._3 < 0))).toMap
      fm.addFactor(new CeFactor(fid, args.map(_._2).toArray, htW(wid), signMap))
    }
    return fm
  }

  // Main parsing function for Ce's format
  def parse(szWeightName: String, szVarName: String, szFactorName: String): FactorManager = {
    val htW = parseWeight(szWeightName)
    return parseVariables(szVarName, szFactorName, htW)
  }
}


// allocate once
// shuffle inside thread
// measure sample time and shuffle time
// larger data sets 100MB