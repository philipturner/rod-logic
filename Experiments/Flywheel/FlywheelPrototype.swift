import Foundation
import HDL
import MM4
import Numerics

// Experiment with a flywheel power source for rod logic.

func createGeometry() -> [[Entity]] {
  let lattice = Lattice<Hexagonal> { h, k, l in
    let h2k = h + 2 * k
    Bounds { 4 * h + 14 * h2k + 20 * l }
    Material { .checkerboard(.carbon, .silicon) }
    
    Volume {
      Origin { 5 * h2k }
      
      Volume {
        Origin { 2 * h + 2 * h2k }
        Concave {
          Plane { -k + h }
          Plane { k + 2 * h }
        }
        Replace { .empty }
      }
      
      Volume {
        Origin { 1 * h2k }
        Plane { -h2k }
        Replace { .atom(.carbon) }
      }
      
      Volume {
        Origin { 3 * h2k }
        Plane { h2k }
        Replace { .atom(.silicon) }
      }
      
      for signY in [Float(-1), 1] {
        Volume {
          Origin { 2 * h2k }
          Origin { signY * 2 * h2k }
          Concave {
            Plane { signY * h2k }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { 18.5 * l }
              Plane { -l }
            }
          }
          Replace { .empty }
          
          for signX in [Float(-1), 1] {
            Concave {
              Plane { signY * h2k }
              if signY == 1 {
                if signX == -1 {
                  Origin { 0.01 * l }
                  Plane { -l }
                } else {
                  Origin { 19.5 * l }
                  Plane { l }
                }
                
                Origin { 1.6 * h2k }
                Plane { -h2k }
                Replace { .atom(.carbon) }
              } else {
                if signX == -1 {
                  
                } else {
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  
  var reconstruction = SurfaceReconstruction()
  reconstruction.topology.insert(atoms: lattice.atoms)
  reconstruction.material = .checkerboard(.carbon, .silicon)
  reconstruction.removePathologicalAtoms()
  reconstruction.createBulkAtomBonds()
  reconstruction.createHydrogenSites()
  reconstruction.resolveCollisions()
  reconstruction.createHydrogenBonds()
  var topology = reconstruction.topology

  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = topology.atoms.map(\.atomicNumber)
  paramsDesc.bonds = topology.bonds
  let parameters = try! MM4Parameters(descriptor: paramsDesc)
  
  var forceFieldDesc = MM4ForceFieldDescriptor()
  forceFieldDesc.parameters = parameters
  let forceField = try! MM4ForceField(descriptor: forceFieldDesc)
  forceField.positions = topology.atoms.map(\.position)
  forceField.minimize()
  
  for i in topology.atoms.indices {
    topology.atoms[i].position = forceField.positions[i]
  }
  
  return [topology.atoms]
  
  var rigidBodyDesc = MM4RigidBodyDescriptor()
  rigidBodyDesc.parameters = parameters
  rigidBodyDesc.positions = topology.atoms.map(\.position)
  var rigidBody1 = try! MM4RigidBody(descriptor: rigidBodyDesc)
  rigidBody1.centerOfMass = .zero
  rigidBody1.rotate(angle: -.pi / 2, axis: [0, 1, 0])
  
  var rigidBodies: [MM4RigidBody] = []
  for rigidBodyID in 0..<3 {
    var rigidBody = rigidBody1
    let angle = Double(rigidBodyID) * 2 * .pi / 3
    rigidBody.rotate(angle: angle, axis: [0, 0, 1])
    
    let quaternion = Quaternion<Double>(angle: angle, axis: [0, 0, 1])
    let displacement = quaternion.act(on: [0, 4.2, 0])
    rigidBody.centerOfMass += displacement
    rigidBodies.append(rigidBody)
  }
  
  var paramsDesc2 = MM4ParametersDescriptor()
  paramsDesc2.atomicNumbers = []
  paramsDesc2.bonds = []
  var parameters2 = try! MM4Parameters(descriptor: paramsDesc2)
  var rigidBodyDesc2 = MM4RigidBodyDescriptor()
  rigidBodyDesc2.positions = []
  for rigidBody in rigidBodies {
    parameters2.append(contentsOf: rigidBody.parameters)
    rigidBodyDesc2.positions! += rigidBody.positions
  }
  rigidBodyDesc2.parameters = parameters2
  var rigidBody2 = try! MM4RigidBody(descriptor: rigidBodyDesc2)
  // n * 1 GHz rotation speed
  rigidBody2.angularMomentum = SIMD3(10 * 0.0062, 0, 0) * rigidBody2.momentOfInertia
  
  
//  return [zip(parameters2.atoms.atomicNumbers, rigidBody2.positions).map {
//    Entity(storage: SIMD4($1, Float($0)))
//  }]
  
  for rigidBodyID in 0..<3 {
    let original = rigidBodies[rigidBodyID]
    var descriptor = MM4RigidBodyDescriptor()
    descriptor.parameters = original.parameters
    descriptor.positions = original.positions
    let span = original.parameters.atoms.count
    let range = (rigidBodyID * span)..<((rigidBodyID + 1) * span)
    descriptor.velocities = Array(rigidBody2.velocities[range])
    rigidBodies[rigidBodyID] = try! MM4RigidBody(descriptor: descriptor)
    print(rigidBodies[rigidBodyID].linearMomentum / rigidBodies[rigidBodyID].mass)
  }
  
  var forceFieldDesc2 = MM4ForceFieldDescriptor()
  forceFieldDesc2.parameters = rigidBody2.parameters
  let forceField2 = try! MM4ForceField(descriptor: forceFieldDesc2)
  forceField2.positions = rigidBodies.flatMap(\.positions)
  
  var animation: [[Entity]] = []
  for frameID in 0...1000 {
    if frameID % 10 == 0 {
      print("frame=\(frameID)")
    }
    if frameID > 0 {
      forceField2.simulate(time: 0.040)
      
      if frameID == 60 || frameID == 70 || frameID == 80 {
        forceField2.velocities = forceField2.velocities.map { $0 * 0 }
      }
      if frameID == 90 {
        forceField2.velocities = rigidBodies.flatMap(\.velocities)
      }
      
//      forceField2.positions = rigidBodies.flatMap(\.positions)
//      forceField2.velocities = rigidBodies.flatMap(\.velocities)
//      let forces = forceField2.forces
//      var cursor = 0
//      let spacing = rigidBodies[0].parameters.atoms.count
//      
//      for rigidBodyID in rigidBodies.indices {
//        let range = cursor..<(cursor + spacing)
//        cursor += spacing
//        var copy = rigidBodies[rigidBodyID]
//        copy.forces = Array(forceField2.forces[range])
//        copy.linearMomentum += 0.040 * copy.netForce!
//        copy.angularMomentum += 0.040 * copy.netTorque!
//        
//        let velocity = copy.linearMomentum / copy.mass
//        let angularVelocity = copy.angularMomentum / copy.momentOfInertia
//        let angularSpeed = (angularVelocity * angularVelocity).sum().squareRoot()
//        copy.centerOfMass += 0.040 * velocity
//        copy.rotate(angle: 0.040 * angularSpeed)
//        rigidBodies[rigidBodyID] = copy
//      }
    }
    
    var frame: [Entity] = []
//    for rigidBodyID in 0..<3 {
//      let rigidBody = rigidBodies[rigidBodyID]
//      for (i, position) in rigidBody.positions.enumerated() {
//        let atomicNumber = rigidBody.parameters.atoms.atomicNumbers[i]
//        let entity = Entity(storage: SIMD4(position, Float(atomicNumber)))
//        frame.append(entity)
//      }
//    }
    for i in rigidBody2.parameters.atoms.indices {
      let atomicNumber = rigidBody2.parameters.atoms.atomicNumbers[i]
      let position = forceField2.positions[i]
      frame.append(Entity(storage: SIMD4(position, Float(atomicNumber))))
    }
    animation.append(frame)
  }
  return animation
}
