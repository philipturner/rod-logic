//
//  GenerateUnit.swift
//  MolecularRendererApp
//
//  Created by Philip Turner on 3/24/24.
//

import Foundation
import HDL
import MM4
import Numerics

struct GenerateUnit {
  // The carry in.
  var carryIn: Rod
  
  // The generate signal.
  //
  // Ordered from bit 0 -> bit 3.
  var signal: [Rod] = []
  
  // The generate signal, transmitted vertically.
  // - keys: The source layer (>= 0), or 'carryIn' (-1).
  // - values: The associated logic rods.
  var probe: [Int: Rod] = [:]
  
  // The carry chains that terminate at the current bit.
  // - keys: The source layer (lane 0) and the destination layer (lane 1).
  // - values: The associated logic rods.
  var broadcast: [SIMD2<Int>: Rod] = [:]
  
  // The rods in the unit, gathered into an array.
  var rods: [Rod] {
    Array([carryIn]) +
    signal +
    Array(probe.values) +
    Array(broadcast.values)
  }
  
  init() {
    // Create 'carryIn'.
    do {
      let z = 5.75 * Float(4)
      let offset = SIMD3(0, 0, z + 0)
      let pattern = GenerateUnit
        .signalPattern(layerID: 0)
      let rod = GenerateUnit
        .createRodX(offset: offset, pattern: pattern)
      carryIn = rod
    }
    
    for layerID in 1...4 {
      let y = 6 * Float(layerID)
      
      // Create 'signal'.
      do {
        let z = 5.75 * Float(4 - layerID)
        let offset = SIMD3(0, y, z + 0)
        let pattern = GenerateUnit
          .signalPattern(layerID: layerID)
        let rod = GenerateUnit
          .createRodX(offset: offset, pattern: pattern)
        signal.append(rod)
      }
      
      // Create 'broadcast'.
      for positionZ in ((4 - layerID) + 1)...4 {
        let z = 5.75 * Float(positionZ)
        let offset = SIMD3(0, y, z + 0)
        let pattern = GenerateUnit
          .broadcastPattern()
        let rod = GenerateUnit
          .createRodX(offset: offset, pattern: pattern)
        
        let key = SIMD2(Int(positionZ), Int(layerID))
        broadcast[key] = rod
      }
    }
    
    // Create 'probe'.
    for positionZ in 0...3 {
      let z = 5.75 * Float(positionZ)
      let offset = SIMD3(10.75, 0, z + 3.25)
      let pattern = GenerateUnit
        .probePattern(positionZ: positionZ)
      let rod = GenerateUnit
        .createRodY(offset: offset, pattern: pattern)
      
      let key = 2 - positionZ
      probe[key] = rod
    }
  }
}

extension GenerateUnit {
  private static func createRodX(
    offset: SIMD3<Float>,
    pattern: KnobPattern
  ) -> Rod {
    let rodLatticeX = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 77 * h + 2 * h2k + 2 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        pattern(h, h2k, l)
      }
    }
    
    let atoms = rodLatticeX.atoms.map {
      var copy = $0
      var position = copy.position
      let latticeConstant = Constant(.square) { .elemental(.carbon) }
      position += SIMD3(0, 0.85, 0.91)
      position += offset * latticeConstant
      copy.position = position
      return copy
    }
    return Rod(atoms: atoms)
  }
  
  private static func createRodY(
    offset: SIMD3<Float>,
    pattern: KnobPattern
  ) -> Rod {
    let rodLatticeY = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 46 * h + 2 * h2k + 2 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        pattern(h, h2k, l)
      }
    }
    
    let atoms = rodLatticeY.atoms.map {
      var copy = $0
      var position = copy.position
      let latticeConstant = Constant(.square) { .elemental(.carbon) }
      position = SIMD3(position.y, position.x, position.z)
      position += SIMD3(0.85, 0, 0.91)
      position += offset * latticeConstant
      copy.position = position
      return copy
    }
    return Rod(atoms: atoms)
  }
}

extension GenerateUnit {
  private static func signalPattern(layerID: Int) -> KnobPattern {
    { h, h2k, l in
      let clockingShift: Float = 4
      
      // Exclude the carry in from these knobs.
      if layerID > 0 {
        // Connect to operand A.
        Volume {
          Concave {
            Convex {
              Origin { 2 * h }
              Plane { h }
            }
            Convex {
              Origin { 0.5 * h2k }
              Plane { -h2k }
            }
            Convex {
              Origin { 7 * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
        
        // Connect to operand B.
        Volume {
          Concave {
            Convex {
              Origin { 11 * h }
              Plane { h }
            }
            Convex {
              Origin { 0.5 * h2k }
              Plane { -h2k }
            }
            Convex {
              Origin { 16 * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
      
      // Exclude the last generate bit from interaction with the probes.
      if layerID != 4 {
        // Create a groove that interacts with the 'probe' rod.
        Volume {
          Concave {
            Convex {
              Origin { (17.5 + clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 0.5 * l }
              Plane { -l }
            }
            Convex {
              Origin { (23 + clockingShift) * h }
              Plane { -h }
            }
          }
          Replace { .empty }
        }
      }
    }
  }
  
  private static func probePattern(positionZ: Int) -> KnobPattern {
    { h, h2k, l in
      let clockingShift: Float = 4
      
      // Groove for interaction with the 'signal' rod.
      do {
        var startOffset: Float
        switch positionZ {
        case 0: startOffset = 28
        case 1: startOffset = 19.5
        case 2: startOffset = 11
        case 3: startOffset = 2.5
        default: fatalError("Unrecognized position Z.")
        }
        let endOffset: Float = 6 + startOffset
        
        Volume {
          Concave {
            Convex {
              Origin { startOffset * h }
              Plane { h }
            }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { endOffset * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
      
      // Groove for interaction with layer 1.
      if positionZ > 2 {
        Volume {
          Concave {
            Convex {
              Origin { (11 - clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { (16.5 - clockingShift) * h }
              Plane { -h }
            }
          }
          Replace { .empty }
        }
      }
      
      // Groove for interaction with layer 2.
      if positionZ > 1 {
        Volume {
          Concave {
            Convex {
              Origin { (19.5 - clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { (25 - clockingShift) * h }
              Plane { -h }
            }
          }
          Replace { .empty }
        }
      }
      
      // Groove for interaction with layer 3.
      if positionZ > 0 {
        Volume {
          Concave {
            Convex {
              Origin { (28 - clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { (33.5 - clockingShift) * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
      
      // Groove for interaction with layer 4.
      do {
        Volume {
          Concave {
            Convex {
              Origin { (36.5 - clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 1.49 * l }
              Plane { l }
            }
            Convex {
              Origin { (42 - clockingShift) * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
    }
  }
  
  // We haven't reached the level of detail where individual broadcasts get
  // unique patterns.
  private static func broadcastPattern() -> KnobPattern {
    { h, h2k, l in
      let clockingShift: Float = 4
      
      // Create a groove to avoid collision with the A/B operands.
      Volume {
        Concave {
          Convex {
            Origin { 0.5 * h2k }
            Plane { -h2k }
          }
          Convex {
            Origin { 16 * h }
            Plane { -h }
          }
          Replace { .empty }
        }
      }
      
      // Create a groove to interact with the 'probe' rod.
      Volume {
        Concave {
          Convex {
            Origin { (21.5 - clockingShift) * h }
            Plane { h }
          }
          Convex {
            Origin { 0.5 * l }
            Plane { -l }
          }
          Convex {
            Origin { (27 - clockingShift) * h }
            Plane { -h }
          }
          Replace { .empty }
        }
      }
    }
  }
}
