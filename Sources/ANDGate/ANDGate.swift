// Save as GitHub Gist instead of polluting the molecular-renderer source tree.

import Foundation
import HDL
import MM4
import Numerics
import OpenMM

// Test two-bit logic gates with full MD simulation. Verify that they work
// reliably at room temperature with the proposed actuation mechanism, at up
// to a 3 nm vdW cutoff. How long do they take to switch?
//
// This may require serializing long MD simulations to the disk for playback.
func createGeometry() -> [Entity] {
  var systemDesc = SystemDescriptor()
  systemDesc.patternA = { h, k, l in
    let h2k = h + 2 * k
    Volume {
      Concave {
        // 2 or 6
        Origin { 6 * h }
        Plane { h }
        
        Origin { 1.5 * h2k }
        Plane { h2k }
        
        // 6 or 6
        Origin { 6 * h }
        Plane { -h }
      }
      
      Replace { .empty }
    }
  }
  systemDesc.patternB = { h, k, l in
    let h2k = h + 2 * k
    Volume {
      Concave {
        // 2 or 6
        Origin { 6 * h }
        Plane { h }
        
        Origin { 1.5 * h2k }
        Plane { h2k }
        
        // 6 or 6
        Origin { 6 * h }
        Plane { -h }
      }
      
      Replace { .empty }
    }
  }
  systemDesc.patternC = { h, k, l in
    let h2k = h + 2 * k
    Volume {
      Concave {
        Origin { 2 * h }
        Plane { h }
        
        Origin { 0.5 * h2k }
        Plane { -h2k }
        
        Origin { 6 * h }
        Plane { -h }
      }
      Concave {
        Origin { 11 * h }
        Plane { h }
        
        Origin { 0.5 * h2k }
        Plane { -h2k }
        
        Origin { 5 * h }
        Plane { -h }
      }
      
      Replace { .empty }
    }
  }
  
  var system = System(descriptor: systemDesc)
  system.passivate()
  system.minimizeSurfaces()
  
  // Initialize the rigid bodies with thermal velocities, then zero out the
  // drift in bulk momentum.
  do {
    let randomThermalVelocities = system.createRandomThermalVelocities()
    system.setAtomVelocities(randomThermalVelocities)
    
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.forces = []
    var rigidBodies = system.createRigidBodies(parametersDescriptor: paramsDesc)
    for rigidBodyID in rigidBodies.indices {
      var rigidBody = rigidBodies[rigidBodyID]
      rigidBody.linearMomentum = .zero
      rigidBody.angularMomentum = .zero
      
      // Give the input drive wall a velocity of 100 m/s in the Y direction.
      //
      // TODO: Rewrite this code. We need to equilibriate the system before
      // giving the drive wall any momentum. In the first iteration, give each
      // atom double the thermal kinetic energy. 40 fs per iteration should be
      // sufficient. Repeat until the energy drift has stabilized.
      if rigidBodyID == 1 {
        let atomVelocitiesBefore = rigidBody.velocities
        let massesBefore = rigidBody.parameters.atoms.masses
        var atomKineticEnergyBefore: Double = .zero
        for atomID in rigidBody.parameters.atoms.indices {
          let mass = massesBefore[atomID]
          let velocityBefore = atomVelocitiesBefore[atomID]
          atomKineticEnergyBefore += Double(0.5 * mass * (velocityBefore * velocityBefore).sum())
        }
        
        let mass = rigidBody.mass
        let velocity = SIMD3<Double>(0, 0.000, 0)
        let linearMomentum = mass * velocity
        rigidBody.linearMomentum = linearMomentum
        
        let kineticEnergy = 0.5 * mass * (velocity * velocity).sum()
        print("Organized kinetic energy:", kineticEnergy)
        
        let atomVelocitiesAfter = rigidBody.velocities
        let massesAfter = rigidBody.parameters.atoms.masses
        var atomKineticEnergyAfter: Double = .zero
        for atomID in rigidBody.parameters.atoms.indices {
          let mass = massesAfter[atomID]
          let velocityAfter = atomVelocitiesAfter[atomID]
          atomKineticEnergyAfter += Double(0.5 * mass * (velocityAfter * velocityAfter).sum())
        }
        
        print("Before:", atomKineticEnergyBefore)
        print("After:", atomKineticEnergyAfter)
      }
      
      rigidBodies[rigidBodyID] = rigidBody
    }
    let rigidBodyAtomVelocities = rigidBodies.map(\.velocities)
    system.setAtomVelocities(rigidBodyAtomVelocities)
  }
  
  let paramsDesc = MM4ParametersDescriptor()
  // use default forces
  let rigidBodies = system.createRigidBodies(parametersDescriptor: paramsDesc)
  
  var emptyParamsDesc = MM4ParametersDescriptor()
  emptyParamsDesc.atomicNumbers = []
  emptyParamsDesc.bonds = []
  var systemParameters = try! MM4Parameters(descriptor: emptyParamsDesc)
  for rigidBodyID in rigidBodies.indices {
    let rigidBody = rigidBodies[rigidBodyID]
    var parameters = rigidBody.parameters
    var boundingBoxMinY: Float = .greatestFiniteMagnitude
    var boundingBoxMaxY: Float = -.greatestFiniteMagnitude
    for atomID in parameters.atoms.indices {
      let position = rigidBody.positions[atomID]
      boundingBoxMinY = min(boundingBoxMinY, position.y)
      boundingBoxMaxY = max(boundingBoxMaxY, position.y)
    }
    if rigidBodyID == 1 {
      boundingBoxMaxY = .signalingNaN
    }
    if rigidBodyID == 2 {
      boundingBoxMinY = .signalingNaN
    }
    
    let frozenRigidBodyIDs: Set<Int> = [0, 2]
    if frozenRigidBodyIDs.contains(rigidBodyID) {
      for atomID in parameters.atoms.indices {
        let centerType = parameters.atoms.centerTypes[atomID]
        guard centerType == .quaternary else {
          continue
        }
        let position = rigidBody.positions[atomID]
        let closenessMin = (position.y - boundingBoxMinY).magnitude
        let closenessMax = (position.y - boundingBoxMaxY).magnitude
        guard closenessMin < 0.5 || closenessMax < 0.5 else {
          continue
        }
        guard Float.random(in: 0..<1) < 0.4 else {
          continue
        }
        parameters.atoms.masses[atomID] = 0
      }
    }
    systemParameters.append(contentsOf: parameters)
  }
  
  let platforms = OpenMM_Platform.platforms
  let reference = platforms.first(where: {
    $0.name == "HIP"
  })!
  
  var forceFieldDesc = MM4ForceFieldDescriptor()
  forceFieldDesc.cutoffDistance = 2.5
  forceFieldDesc.integrator = .multipleTimeStep
  forceFieldDesc.parameters = systemParameters
  forceFieldDesc.platform = reference
  let forceField = try! MM4ForceField(descriptor: forceFieldDesc)
  let systemStartPositions = rigidBodies.flatMap(\.positions)
  let systemStartVelocities = rigidBodies.flatMap(\.velocities)
  forceField.positions = systemStartPositions
  
  do {
    let kinetic = forceField.energy.kinetic
    let potential = forceField.energy.potential
    print("frame:", -1, terminator: " | ")
    print("kinetic:", kinetic, terminator: " | ")
    print("potential:", potential, terminator: " | ")
  }
  
  forceField.minimize()
  forceField.velocities = systemStartVelocities
  
  var frames: [[Entity]] = []
  var initialEnergy: Double?
  for frameID in 0...2 {
    if frameID > 0 {
      forceField.simulate(time: 0.040)
    }
    let kinetic = forceField.energy.kinetic
    let potential = forceField.energy.potential
    print("frame:", frameID, terminator: " | ")
    print("kinetic:", kinetic, terminator: " | ")
    print("potential:", potential, terminator: " | ")
    if frameID > 0 {
      let totalEnergy = kinetic + potential
      let driftEnergy = totalEnergy - initialEnergy!
      print("drift:", driftEnergy)
    } else {
      initialEnergy = kinetic + potential
    }
    
    var frame: [Entity] = []
    for atomID in systemParameters.atoms.indices {
      let atomicNumber = systemParameters.atoms.atomicNumbers[atomID]
      let position = forceField.positions[atomID]
      let element = Element(rawValue: atomicNumber)!
      let entity = Entity(position: position, type: .atom(element))
      frame.append(entity)
    }
    frames.append(frame)
  }
  
  // Next, give the drive walls velocities. Send them through the clocking
  // motion at 100 m/s and time the gate switching time.
  
  // TODO: Save the above code somewhere. Do an RBD simulation before doing
  // the MD simulation.
  //
  // Nevermind. Rigid body dynamics is easier to set up in general. I'll
  // deal with the issues with setting up MD simulations another time.
  
  return frames
}
