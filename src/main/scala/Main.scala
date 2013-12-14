package org.hazy.gibbs

import java.io.File

object GibbsSampler extends App {

  // Parse command line arguments
  case class Config(variablesFile: File, factorsFile: File, weightsFile: File, 
    numSamples: Int, numThreads: Int, outputFile: File)

  val parser = new scopt.OptionParser[Config]("gibbs_sampler") {
    head("gibbs_sampler", "")
    opt[File]('v', "variables") required() valueName("<file>") action { (x, c) =>
      c.copy(variablesFile = x) } text("variables file (required)")
    opt[File]('f', "factors") required() valueName("<file>") action { (x, c) =>
      c.copy(factorsFile = x) } text("factors file (required) ")
    opt[File]('w', "weights") required() valueName("<file>") action { (x, c) =>
      c.copy(weightsFile = x) } text("factors file (required)")
    opt[File]('o', "output") required() valueName("<file>") action { (x, c) =>
      c.copy(outputFile = x) } text("output file (required)")
    opt[Int]('n', "samples") required() action { (x, c) =>
      c.copy(numSamples = x) } text("number of samples (required)")
    opt[Int]('t', "threads") required() action { (x, c) =>
      c.copy(numThreads = x) } text("number of threads (required)")
  }

  val initialConfig = Config(null, null, null, 100, 0, null)
  parser.parse(args, initialConfig) map { config =>
    val fm = DBParser.parse(config.weightsFile.getCanonicalPath, 
      config.variablesFile.getCanonicalPath,
      config.factorsFile.getCanonicalPath)
    fm.MarginalsFile(config.numSamples, fm.variables, config.numThreads, 
      config.outputFile.getCanonicalPath)
  } getOrElse {
    System.exit(1)
  }

}