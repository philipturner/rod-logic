// Save as GitHub Gist instead of polluting the molecular-renderer source tree.

import Foundation
import HDL
import MM4
import Numerics

// Create a piece of code that detects collisions during a scripted animation.
// Find which structures are colliding, and highlight them in red. Log the pairs
// of structures that have collided.
func createGeometry() -> [[Entity]] {
  // Create the initial blocks.
  let originalBlock = Block()
  
  var block1 = originalBlock
  block1.rigidBody.centerOfMass += SIMD3(-2, 0, 0)
  
  var block2 = originalBlock
  var block2Axis = SIMD3<Double>(1, 2, 3)
  block2Axis /= (block2Axis * block2Axis).sum().squareRoot()
  block2.rigidBody.rotate(angle: 1, axis: block2Axis)
  block2.rigidBody.centerOfMass += SIMD3(0, 0.2, 0)
  
  var block3 = block2
  block3.rigidBody.centerOfMass += SIMD3(2, 0, 0)
  
  var block4 = block3
  block4.rigidBody.centerOfMass += SIMD3(2, 0, 0)
  
  // Loop over the frames in an animation.
  var frames: [[Entity]] = []
  for frameID in 0..<600 {
    // Update moving blocks this frame.
    block1.rigidBody.centerOfMass += SIMD3(0.02, 0, 0)
    
    // Collect the blocks into an array.
    var blocks = [block1, block2, block3, block4]
    for blockID in blocks.indices {
      var block = blocks[blockID]
      for i in block.topology.atoms.indices {
        let position = block.rigidBody.positions[i]
        block.topology.atoms[i].position = position
      }
      blocks[blockID] = block
    }
    
    // Transform the atoms into a list, check for per-atom collisions globally.
    // Only check for whether carbon atoms in different blocks are close.
    var carbons: [Entity] = []
    var carbonsToBlocksMap: [Int] = []
    for blockID in blocks.indices {
      let block = blocks[blockID]
      for atom in block.topology.atoms {
        if atom.atomicNumber == 6 {
          carbons.append(atom)
          carbonsToBlocksMap.append(blockID)
        }
      }
    }
    
    // Mark colliding blocks as red, and report pairs of overlapping blocks.
    var collisionPairs: [SIMD2<Int>: Bool] = [:]
    var collisionTopology = Topology()
    collisionTopology.insert(atoms: carbons)
    let matches = collisionTopology.match(
      carbons, algorithm: .covalentBondLength(1.5), maximumNeighborCount: 16)
    for lhsAtomID in carbons.indices {
      for rhsAtomID in matches[lhsAtomID] {
        let lhsBlockID = carbonsToBlocksMap[lhsAtomID]
        let rhsBlockID = carbonsToBlocksMap[Int(rhsAtomID)]
        guard lhsBlockID != rhsBlockID else {
          continue
        }
        let pair = SIMD2(
          min(lhsBlockID, rhsBlockID),
          max(lhsBlockID, rhsBlockID))
        collisionPairs[pair] = true
      }
    }
    var redBlocksMask = [Bool](repeating: false, count: blocks.count)
    for pair in collisionPairs.keys {
      print("frames[\(frameID)] - collision (\(pair.x), \(pair.y))")
      redBlocksMask[pair.x] = true
      redBlocksMask[pair.y] = true
    }
    for blockID in blocks.indices {
      if redBlocksMask[blockID] {
        var block = blocks[blockID]
        for atomID in block.topology.atoms.indices {
          block.topology.atoms[atomID].atomicNumber = 8
        }
        blocks[blockID] = block
      }
    }
    
    // Append the frame.
    let frame = blocks.flatMap { $0.topology.atoms }
    frames.append(frame)
  }
  
  return frames
}

// Compile a rectangular block, for testing a collision detection engine.
struct Block {
  var topology: Topology
  var rigidBody: MM4RigidBody!
  
  // This initializer has a very high amount of latency. Call it sparingly.
  init() {
    topology = Topology()
    createLattice()
    reconstructSurface()
    minimizeSurface()
    createRigidBody()
  }
  
  mutating func createLattice() {
    let blockLattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 4 * h + 2 * h2k + 2 * l }
      Material { .elemental(.carbon) }
    }
    topology.insert(atoms: blockLattice.atoms)
  }
  
  mutating func reconstructSurface() {
    var reconstruction = SurfaceReconstruction()
    reconstruction.material = .elemental(.carbon)
    reconstruction.topology = topology
    reconstruction.removePathologicalAtoms()
    reconstruction.createBulkAtomBonds()
    reconstruction.createHydrogenSites()
    reconstruction.resolveCollisions()
    reconstruction.createHydrogenBonds()
    topology = reconstruction.topology
    topology.sort()
  }
  
  mutating func minimizeSurface() {
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = topology.atoms.map(\.atomicNumber)
    paramsDesc.bonds = topology.bonds
    var parameters = try! MM4Parameters(descriptor: paramsDesc)
    for i in topology.atoms.indices {
      if parameters.atoms.centerTypes[i] == .quaternary {
        parameters.atoms.masses[i] = 0
      }
    }
    
    var forceFieldDesc = MM4ForceFieldDescriptor()
    forceFieldDesc.parameters = parameters
    let forceField = try! MM4ForceField(descriptor: forceFieldDesc)
    forceField.positions = topology.atoms.map(\.position)
    forceField.minimize()
    
    for i in topology.atoms.indices {
      topology.atoms[i].position = forceField.positions[i]
    }
  }
  
  mutating func createRigidBody() {
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = topology.atoms.map(\.atomicNumber)
    paramsDesc.bonds = topology.bonds
    paramsDesc.forces = [.nonbonded]
    let parameters = try! MM4Parameters(descriptor: paramsDesc)
    
    var rigidBodyDesc = MM4RigidBodyDescriptor()
    rigidBodyDesc.parameters = parameters
    rigidBodyDesc.positions = topology.atoms.map(\.position)
    rigidBody = try! MM4RigidBody(descriptor: rigidBodyDesc)
  }
}
