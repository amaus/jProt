# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Section on Semantic Versioning to README.
- Comments in build.gradle to explain Semantic Versioning.
- Comments in build.gradle to explain gradle build and install usage.
- class AminoAcid
- class Bond
- class EncadParameters
- class Molecule
### Changed
- Internal Change: changed contains(String flag) in CommandLineParser to no
  longer split on colons. The proper format for arguments if `-flag value`
### Deprecated
[None]
### Removed
[None]
### Fixed
[None]
### Security
[None]

## [0.5.0] - 2017-7-17
### Added
- Section header in README for the blurb about the jMath Dependency
- Releases all the way through [0.5.0] to this CHANGELOG
### Changed
- Renamed class MausMetrics to Metrics and executable AngularDistanceJProt to
  JProtMetrics.
### Removed
- Executable PointRunner from the project.

## [0.4.4] - 2017-6-7
### Changed
- The Regions of Local Similarity now have a default threshold of 1.0 Angstroms.
  The percent of residues in the top 4 regions is now printed out.
### Fixed
- fixed type in AngularDistanceJProt when setting the GDT thresholds.

## [0.4.3] - 2017-5-11
### Changed
- updated Global Distance Test Thresholds in AngularDistanceJProt from 2, 4,
  and 8 to 1, 2, 4, and 8 with high accuracy threshold at half of those.

## [0.4.2] - 2017-4-25
### Changed
- updated MausMetrics to use a faster MAX CLIQUE solver from jMath. The
  speed up is significant.

## [0.4.1] - 2017-2-23
### Changed
- updated logic in MausMetrics when finding clique covering. Removed the clique
  covering logic because it has been added to jMath.
- minor output fix in AngularDistanceJProt

## [0.4.0] - 2016-11-20
### Changed
- AngularDistanceJProt no longer accepts a distance differences file as input.
  As before, it requires a distances file for each structure. This is a
  backwards incompatible change, but since this project is not at the production
  stage, I am not updating the MAJOR version number. See semver.org for details.
### Fixed
- Fixed bug in MausMetrics when printing out pymol script. Now prints out the
  proper residue IDs for the respective structures in the pymol selection
  statements.

## [0.3.0] - 2016-11-17
### Added
- added class CommandLineParser to help with interpreting command line arguments
  in executables. Included a constructor that takes the String[] args,
  a contains(String flag) method, and a getValue(String flag) method.
### Changed
- Big changes in AngularDistanceJProt. It now takes dynamic cmd line args to
  specify which metrics and thresholds to use. Will also print out a lengthy
  help message when called with no arguments or with the -h flag.
- Updated MausMetrics to be more robust when dealing with distance matrices.
- minor updates to javadoc comments
### Fixed
- fixed javadoc compiler errors on method parameter names in various classes.

## [0.2.0] - 2016-11-12
### Added
- new class MausMetrics with AngularDistance, Covering Regions of Local
  Similarity, and Global Distance Test using Max Cliques (percent of residues
  in a clique built using distance threshold of 1.0 Angstroms, then percent
  of residues in clique for 2.0 A, and 4.0 A)
- Added Executable AngularDistanceJProt to run these metrics and print pymol
  scripts to color residues in these regions
- Added explanation in README about the jMath dependency
### Changed
- Project name to jProt instead of jProtein
- cleaned up build.gradle
- updated .gitignore to ignore gradle.properties
### Deprecated
- Executable PointRunner

## [0.1.0] - 2016-11-11
### Added
- initial build.gradle file with dependencies, the ability to create executables
  scripts, and basic information on project: group, version, plugins (maven,
  maven-publish, application)
- README.md
- .gitignore
- Atom class with constructors, getCoordinates(), distance(Atom other),
  moveTo(double x, double y, double z), and toString()
- PointRunner executable, first basic executable for testing purposes. Will not
  remain in project.
