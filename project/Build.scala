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


import java.io.File
import sbt._
import Keys._
import sbt.classpath.ClasspathUtilities


object GenomeWeaverBuild extends Build {


  val SCALA_VERSION = "2.9.2"

  lazy val buildSettings = Defaults.defaultSettings ++ Seq[Setting[_]](
    organization := "org.xerial",
    organizationName := "Xerial Project",
    organizationHomepage := Some(new URL("http://xerial.org/")),
    version := "0.1-SNAPSHOT",
    description := "Genome Weaver: Toolkit for Genome Sciences",
    scalaVersion := SCALA_VERSION,
    publishMavenStyle := true,
    publishArtifact in Test := false,
    publishTo <<= version {
      v: String =>
        val repoPath = "/home/web/maven.xerial.org/repository/" + (if (v.trim.endsWith("SNAPSHOT")) "snapshot" else "artifact")
        Some(Resolver.ssh("Xerial Repo", "www.xerial.org", repoPath)
          as(System.getProperty("user.name"), new File(Path.userHome.absolutePath + "/.ssh/id_rsa")))
    },
    otherResolvers := Seq(Resolver.file("localM2", file(Path.userHome.absolutePath + "/.m2/repository"))),
    publishLocalConfiguration <<= (packagedArtifacts, deliverLocal, checksums, ivyLoggingLevel) map {
      (arts, _, cs, level) => new PublishConfiguration(None, "localM2", arts, cs, level)
    },
    pomIncludeRepository := {
      _ => false
    },
    parallelExecution := true,
    crossPaths := false,
    resolvers ++= Seq("Typesafe repository" at "http://repo.typesafe.com/typesafe/releases",
      "UTGB Maven repository" at "http://maven.utgenome.org/repository/artifact/",
      "Xerial Maven repository" at "http://www.xerial.org/maven/repository/artifact",
      "Local Maven repository" at "file://" + Path.userHome.absolutePath + "/.m2/repository"
    ),
    scalacOptions ++= Seq("-encoding", "UTF-8", "-deprecation", "-unchecked"),
    pomExtra := {
      <licenses>
        <license>
          <name>Apache 2</name>
          <url>http://www.apache.org/licenses/LICENSE-2.0.txt</url>
        </license>
      </licenses>
        <scm>
          <connection>scm:git:github.com/xerial/genome-weaver.git</connection>
          <developerConnection>scm:git:git@github.com:xerial/genome-weaver.git</developerConnection>
        </scm>
        <properties>
          <scala.version>2.9.2</scala.version>
          <project.build.sourceEncoding>UTF-8</project.build.sourceEncoding>
        </properties>
        <developers>
          <developer>
            <id>leo</id>
            <name>Taro L. Saito</name>
            <url>http://xerial.org/leo</url>
          </developer>
        </developers>
    }
  )

  val distAllClasspaths = TaskKey[Seq[Classpath]]("dist-all-classpaths")
  val distDependencies = TaskKey[Seq[File]]("dist-dependencies")
  val distLibJars = TaskKey[Seq[File]]("dist-lib-jars")

  lazy val distSettings: Seq[Setting[_]] = Seq(
    distAllClasspaths <<= (thisProjectRef, buildStructure) flatMap allProjects(dependencyClasspath.task in Compile),
    distDependencies <<= distAllClasspaths map {
      _.flatten.map(_.data).filter(ClasspathUtilities.isArchive).distinct
    },
    distLibJars <<= (thisProjectRef, buildStructure) flatMap allProjects(packageBin.task in Compile)
  )

  def allProjects[T](task: SettingKey[Task[T]])(currentProject: ProjectRef, structure: Load.BuildStructure): Task[Seq[T]] = {
    val projects: Seq[String] = currentProject.project +: childProjectNames(currentProject, structure)
    projects flatMap {
      task in LocalProject(_) get structure.data
    } join
  }

  def childProjectNames(currentProject: ProjectRef, structure: Load.BuildStructure): Seq[String] = {
    val children = Project.getProject(currentProject, structure).toSeq.flatMap(_.aggregate)
    children flatMap {
      ref =>
        ref.project +: childProjectNames(ref, structure)
    }
  }

  object Dependencies {

    val classWorld = "org.codehaus.plexus" % "plexus-classworlds" % "2.4"

    val testLib = Seq(
      "junit" % "junit" % "4.10" % "test",
      "org.scalatest" %% "scalatest" % "2.0.M1" % "test"
      //"org.hamcrest" % "hamcrest-core" % "1.3.RC2" % "test"
    )

    val bootLib = Seq(
      classWorld
    )

    val coreLib = Seq(
      "org.scala-lang" % "scalap" % SCALA_VERSION,
      "org.xerial" % "xerial-core" % "2.1",
      "org.utgenome" % "utgb-core" % "1.5.8",
      //"org.xerial.silk" % "silk-core" % "0.4",
      "org.javassist" % "javassist" % "3.15.0-GA"
      //"io.netty" % "netty" % "3.3.0.Final"
      //"com.typesafe.akka" % "akka-actor" % "2.0",
      //"com.typesafe.akka" % "akka-remote" % "2.0"
    )
  }


  private val dependentScope = "test->test;compile->compile"

  import Dependencies._

  lazy val root = Project(
    id = "genome-weaver",
    base = file("."),
    settings = buildSettings ++ distSettings ++ Release.settings
      ++ Seq(packageDistTask)
      ++ Seq(libraryDependencies ++= bootLib ++ testLib ++ coreLib)
  ) aggregate(silk, gwLens, gwAlign) dependsOn (silk)

  lazy val silk = RootProject(file("silk"))

  lazy val gwLens = Project(
    id = "lens",
    base = file("lens"),
    settings = buildSettings
      ++ Seq(libraryDependencies +=
      "org.apache.commons" % "commons-compress" % "1.4.1"
    )
  ) dependsOn (silk % dependentScope)


  lazy val gwAlign = Project(
    id = "align",
    base = file("align"),
    settings = buildSettings
      ++ Seq(libraryDependencies ++= bootLib ++ testLib ++ coreLib)
  ) dependsOn (gwLens % dependentScope)


  lazy val copyDependencies = TaskKey[Unit]("copy-dependencies")

  def copyDepTask = copyDependencies <<= (update, crossTarget, scalaVersion) map {
    (updateReport, out, scalaVer) =>
      updateReport.allFiles foreach {
        srcPath =>
          val destPath = out / "lib" / srcPath.getName
          IO.copyFile(srcPath, destPath, preserveLastModified = true)
      }
  }

  lazy val packageDist: TaskKey[File] = TaskKey[File]("package-dist")

  def packageDistTask: Setting[Task[File]] = packageDist <<= (update, version, distLibJars, distDependencies, streams, target, dependencyClasspath in Runtime, classDirectory in Compile, baseDirectory) map {
    (up, ver, libs, depJars, out, target, dependencies, classDirectory, base) => {

      val distDir = target / "dist"

      out.log.info("output dir: " + distDir)
      IO.delete(distDir)
      distDir.mkdirs()

      out.log.info("Copy libraries")
      val libDir = distDir / "lib"
      libDir.mkdirs()
      (libs ++ depJars).foreach(l => IO.copyFile(l, libDir / l.getName))

      out.log.info("Create bin folder")
      val binDir = distDir / "bin"
      binDir.mkdirs()
      IO.copyDirectory(base / "src/script", binDir)

      out.log.info("Generating version info")
      IO.write(distDir / "VERSION", ver)
      out.log.info("done.")

      distDir
    }
  }


}








