# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## Unreleased
### Added
[None]
### Changed
[None]
### Deprecated
[None]
### Removed
[None]
### Fixed
[None]
### Security
[None]

## [0.7.0] - 2018-02-26
### Added
- package sequence
  -  Class Sequence w/ subclasses ProteinSequence, DNASequence, RNASequence
    - sequences of subclass types can be aligned
  - Class Alignment
    - Holds the results of aligning two sequences
- Executable AlignSequences that can align the sequences of 2 PDBs or FASTA
  files.
- Constructor to Residue that takes a char to indicate which amino acid to build
- Feature to PDDFileIO: can now write out bare bones PDB files.
- Residue of Type UNK to residue class and templates.
  - 3 letter name: UNK, 1 letter name: X, full name Unknown
  - Template Coordinates and bonds: same as Alanine
- applyTransformation() methods to Protein, PolypeptideChain, Residue, and Atom
  - applyTransformation() naively applies the transformation to all components
    of that class. Eg, when called on Protein, it calls the applyTransformation
    for each PolypeptideChain, which does the same for every Residue, which then
    applies the transformation to each Atom - whose coordinates are then
    transformed.
- getOmegaAngle(), getPhiAngle(), and getPsiAngle() in PolypeptideChain. Each
  takes a residueID and returns the angle if it can be calculated or 1000
  otherwise.
  - a test to TestProtein to check that all the angles for 1rop are calculated
    correctly.
- setOmegaAngle(), setPhiAngle, and setPsiAngle() in PolypeptideChain
- In Residue, the ability to get and set any of the dihedral angles, three-atom
  angles, and bond lengths.
- Initial implementation for VirtualRibosome, constructs a single chain with
  correct bond angles and residues on opposite sides of the backbone.
- JProtMetrics: new cmd line arguments
  - --chimera to specify to print out chimera script instead of pymol script.
  - --gdt-plot to and prints out the percent of residues for each threshold in
    the range 0.5 to 10.0 in 0.5 increments. This data can be used to produce
    GDT Plots like those on prediction center.
### Changed
- Moved SequenceAligner into new package sequence
- Protein and PolypeptideChain getSequence methods to use new sequence package
- All structure classes are now transformable.
- Atom constructors now take doubles for the coordinates rather than Strings.
- In Metrics, changed the score calculation to be based on the number of
  residues in the reference structure rather than the number of Residues
  in the alignment. This makes the scores more stable and makes more
  sense when evaluating models. If a model only includes a fraction of the
  residues that are in the reference, it should be penalized for the
  missing residues.
### Deprecated
- JProtMetrics no longer accepts csv files as input. Two pdb files must be
  provided.
### Removed
- method getBondLength() and getBondStrength() in Bond. replaced with
  getLength() and getStrength().
- all methods to do with enabling or disabling hydrogens in a protein. This was
  too much extra complication. A protein by default will have all the atoms
  included in the pdb if read in or all heavy atoms and hydrogens is built
  synthetically.

## [0.6.4] 2017-10-27
### Fixed
- PDBIO atom charges bug. Reading in the charges is more flexible. Previously
  only accepted input in form specified by PDB File Format V. 3.30. Now reads in
  charges in form 2+, 1-, +2, -1, 3, etc.

## [0.6.3] 2017-10-16
### Fixed
- Issue-#003: Fix Bond::hashCode(). Mod the final number by MAXINT to keep it in
  range.

## [0.6.2] 2017-10-02
### Added
- A check when adding disulfide bridges to ensure that both chains and residues
  are present in the protein before attempting to add the disulfide bridge.
### Fixed
- BUG issue-#002: Calculating CA Distance Matrices when residues are missing CA.
  A residue is only constructed and added to a protein if the N, CA, and C atoms
  are present.

## [0.6.1] 2017-10-01
### Fixed
- BUG issue-#001: Reading in PDBs without explicit element fields in the Atom
  records. The atomName (CA, CG1, N, O, etc) is now used to determine the
  element.

## [0.6.0] - 2017-10-01
### Added
- Section on Semantic Versioning to README.
- Comments in build.gradle to explain Semantic Versioning.
- Comments in build.gradle to explain gradle build and install usage.
- class Residue
  - JUnit test class for Residue. It tests the construction of several default
    residues.
- Residue data files in main resources. These files specify atom bonds and
  alternate atom type IDs for various energy functions.
- class Bond
- class EncadParameters
- class Molecule
- class PolypeptideChain
- class Protein
- class PDBFileIO with initial ability to read in PDBs. It does not fill in
  missing atoms, but it can build a protein with all chains, residues, and
  atoms specified in the file. It can not write proteins back out to PDB yet.
- stubbed test classes TestProtein TestPolypeptideChain along with
  resource file 1rop.pdb
- SequenceAligner class that can align Protein, DNA, or RNA sequences using the
  Needleman-Wuncsh algorithm with Affine Gap Penalties.
  - Resource files for the Match score matrices, BLOSUM62, DNA, and RNA. Other
    BLOSUM matrices (and PAM) can be added in the future. The DNA and RNA
    matrices are naive and can certainly be improved.
  - junit tests for SequenceAligner
- Add PMF Classes
  - PMFFileIO to handle reading and writing from the encad PMF format
  - PotentialOfMeanForce to represent a PMF with getters and setters for
    its various information. Includes methods to get and set the
    energy functions.
- DistanceMatrixCalculator class in tools package with JUnit tests for base
  and masked versions of the calculateDistanceMatrix methods.
  - This class contains the methods that were formally in Protein and
    PolypeptideChain. The code was duplicate so it has been moved into its
    own standalone static class.
### Changed
- Internal Change: changed contains(String flag) in CommandLineParser to no
  longer split on colons. The proper format for arguments if `-flag value`
- JProtMetrics can now read in two PDBs to specify the proteins to be compared.
  It will perform a sequence alignment and perform the comparison on the
  residues that were aligned.
- Move the usage info for JProtMetrics into a resources file. JProtMetrics reads
  that file in to display usage. All executable will handle usage info this way.
  This design makes editing the usage and handling formatting much easier.
### Deprecated
- JProtMetrics usage ability to take in CA Distance Matrices. In the future,
  it will only be able to work with PDBs as input.
### Removed
- Atom::moveTo()
- Bond:getEnergy()
### Fixed
- Dates in CHANGELOG to better conform to ISO 8601
- Javadoc generation to include links to standard java APIs
### Security
[None]

## [0.5.0] - 2017-07-17
### Added
- Section header in README for the blurb about the jMath Dependency
- Releases all the way through [0.5.0] to this CHANGELOG
### Changed
- Renamed class MausMetrics to Metrics and executable AngularDistanceJProt to
  JProtMetrics.
### Removed
- Executable PointRunner from the project.

## [0.4.4] - 2017-06-07
### Changed
- The Regions of Local Similarity now have a default threshold of 1.0 Angstroms.
  The percent of residues in the top 4 regions is now printed out.
### Fixed
- fixed type in AngularDistanceJProt when setting the GDT thresholds.

## [0.4.3] - 2017-5-11
### Changed
- updated Global Distance Test Thresholds in AngularDistanceJProt from 2, 4,
  and 8 to 1, 2, 4, and 8 with high accuracy threshold at half of those.

## [0.4.2] - 2017-04-25
### Changed
- updated MausMetrics to use a faster MAX CLIQUE solver from jMath. The
  speed up is significant.

## [0.4.1] - 2017-02-23
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
