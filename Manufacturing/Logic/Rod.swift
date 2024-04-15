//
//  Rod.swift
//  MolecularRendererApp
//
//  Created by Philip Turner on 4/11/24.
//

import Foundation
import HDL
import MM4
import Numerics

struct Rod {
  var rigidBody: MM4RigidBody
  
  init(lattice: Lattice<Hexagonal>) {
    let topology = Self.createTopology(lattice: lattice)
    rigidBody = Self.createRigidBody(topology: topology)
  }
  
  static func createTopology(lattice: Lattice<Hexagonal>) -> Topology {
    var reconstruction = SurfaceReconstruction()
    reconstruction.material = .elemental(.carbon)
    reconstruction.topology.insert(atoms: lattice.atoms)
    reconstruction.compile()
    
    var topology = reconstruction.topology
    let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
    
    var removedAtoms: [UInt32] = []
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      guard atom.atomicNumber == 1 else {
        continue
      }
      for j in atomsToAtomsMap[i] {
        let other = topology.atoms[Int(j)]
        if other.atomicNumber == 15 {
          removedAtoms.append(UInt32(i))
        }
        if other.atomicNumber == 16 {
          removedAtoms.append(UInt32(i))
        }
      }
    }
    topology.remove(atoms: removedAtoms)
    topology.sort()
    return topology
  }
  
  static func createRigidBody(topology: Topology) -> MM4RigidBody {
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = topology.atoms.map(\.atomicNumber)
    paramsDesc.bonds = topology.bonds
    let parameters = try! MM4Parameters(descriptor: paramsDesc)
    
    var rigidBodyDesc = MM4RigidBodyDescriptor()
    rigidBodyDesc.parameters = parameters
    rigidBodyDesc.positions = topology.atoms.map(\.position)
    return try! MM4RigidBody(descriptor: rigidBodyDesc)
  }
}
