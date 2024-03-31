import Foundation
import HDL
import MM4
import Numerics
import OpenMM

// WARNING: The renderer could be in 'MRSceneSize.extreme'. If so, it will not
// render any animations.
func createGeometry() -> [[Entity]] {
  // Create the housing.
  let testHousing = TestHousing()
  
  // Instantiate the circuit, then select to part to examine.
  let circuit = Circuit()
  let rod = circuit.propagate.signal[0]
  var testRod = TestRod(rod: rod)
  testRod.rigidBody.rotate(angle: .pi / 2, axis: SIMD3(0, 1, 0))
  testRod.rigidBody.centerOfMass.z += 6
  
  // Instantiate the test system.
  var testSystem = TestSystem(testHousing: testHousing, testRod: testRod)
  testSystem.minimize()
  testSystem.equilibriate(temperature: 298)
  
  // Conduct an ambient MD simulation.
  var frames: [[Entity]] = []
  frames.append(testSystem.createFrame())
  print("frame: 0")
  for frameID in 1...1000 {
    // Update the force field's positions.
    var positions: [SIMD3<Float>] = []
    positions += testSystem.testHousing.rigidBody.positions
    positions += testSystem.testRod.rigidBody.positions
    testSystem.forceField.positions = positions
    
    // Update the force field's velocities.
    var velocities: [SIMD3<Float>] = []
    velocities += testSystem.testHousing.rigidBody.velocities
    velocities += testSystem.testRod.rigidBody.velocities
    testSystem.forceField.velocities = velocities
    
    // Simulate for one timestep.
    testSystem.forceField.simulate(time: 0.040)
    testSystem.updateRigidBodies()
    
    // Record the frame for rendering.
    frames.append(testSystem.createFrame())
    print("frame:", frameID, terminator: " | ")
    
    // Report the time, assuming the timestep is 0.040 fs.
    let time = Double(frameID) * 0.040
    let timeRepr = String(format: "%.1f", time)
    print("time =", timeRepr, "ps", terminator: " | ")
    
    let kineticEnergy = testSystem.forceField.energy.kinetic
    let kineticEnergyRepr = String(format: "%.1f", kineticEnergy)
    print("K =", kineticEnergyRepr, "zJ", terminator: " | ")
    
    // Report how far the pieces have moved.
    let housingRigidBody = testSystem.testHousing.rigidBody!
    let rodRigidBody = testSystem.testRod.rigidBody
    var delta = SIMD3<Float>(housingRigidBody.centerOfMass)
    delta -= SIMD3<Float>(rodRigidBody.centerOfMass)
    let distance = (delta * delta).sum().squareRoot()
    let distanceRepr = String(format: "%.3f", distance)
    print("separation =", distanceRepr, "nm")
  }
  return frames
}