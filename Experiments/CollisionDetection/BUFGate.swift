// Save as GitHub Gist instead of polluting the molecular-renderer source tree.

import Foundation
import HDL
import MM4
import Numerics

// Create a BUF/NOT gate simulator, where collision detection is used to
// advance the output rod. Skip surface passivation and energy minimization.
func createGeometry() -> [[Entity]] {
  let lattice = Lattice<Hexagonal> { h, k, l in
    let h2k = h + 2 * k
    Bounds { 15 * h + 2 * h2k + 2 * l }
    Material { .elemental(.carbon) }
    
    Volume {
      Origin { 1.5 * h2k }
      Concave {
        Plane { h2k }
        Origin { 3 * h }
        Plane { -h }
      }
      Concave {
        Plane { h2k }
        Origin { 7 * h }
        Plane { h }
      }
      Replace { .empty }
    }
  }
  
  // Set the starting position for the input rod.
  var atomsInput = lattice.atoms
  for atomID in atomsInput.indices {
    var atom = atomsInput[atomID]
    atom.position.y = -atom.position.y
    atom.position.y += 2.0
    atom.position = SIMD3(atom.position.z, atom.position.y, atom.position.x)
    atom.position.z += -0.9
    atomsInput[atomID] = atom
  }
  
  #if false
  // Change into a NOT gate.
  var inputCenterOfMass: SIMD3<Float> = .zero
  for atom in atomsInput {
    inputCenterOfMass += atom.position
  }
  inputCenterOfMass /= Float(atomsInput.count)
  for atomID in atomsInput.indices {
    var atom = atomsInput[atomID]
    var deltaZ = atom.position.z - inputCenterOfMass.z
    deltaZ = -deltaZ
    atom.position.z = deltaZ + inputCenterOfMass.z
    atomsInput[atomID] = atom
  }
  #endif
  
  // Set the starting position for the output rod.
  var atomsOutput = lattice.atoms
  for atomID in atomsOutput.indices {
    var atom = atomsOutput[atomID]
    atom.position.x += -2.3
    atomsOutput[atomID] = atom
  }
  
  var frames: [[Entity]] = []
  for frameID in 0..<600 {
    var nextAtomsOutput = atomsOutput
    for atomID in nextAtomsOutput.indices {
      nextAtomsOutput[atomID].position.x += 0.02
    }
    
    // The list of matches corresponds to atoms passed as the array parameter
    // to 'match'.
    var inputTopology = Topology()
    inputTopology.insert(atoms: atomsInput)
    let matches = inputTopology.match(
      atomsOutput, algorithm: .covalentBondLength(3))
    
    var foundMatch = false
    for outputAtomID in atomsOutput.indices {
      if matches[outputAtomID].count > 0 {
        foundMatch = true
      }
    }
    if foundMatch {
      // Don't advance the output rod.
    } else {
      atomsOutput = nextAtomsOutput
    }
    
    frames.append(atomsInput + atomsOutput)
  }
  
  return frames
}
