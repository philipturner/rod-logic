//
//  PropagateUnit.swift
//  MolecularRendererApp
//
//  Created by Philip Turner on 3/24/24.
//

import Foundation
import HDL
import MM4
import Numerics

struct PropagateUnit {
  // The propagate signal.
  //
  // Ordered from bit 0 -> bit 3.
  var signal: [Rod] = []
  
  // The propagate signal, transmitted vertically.
  // - keys: The source layer.
  // - values: The associated logic rods.
  var probe: [Int: Rod] = [:]
  
  // The propagate signal, broadcasted to every applicable carry chain.
  // - keys: The source x-index (0) and the destination layer (1).
  // - values: The associated logic rods.
  var broadcast: [SIMD2<Int>: Rod] = [:]
  
  // The rods in the unit, gathered into an array.
  var rods: [Rod] {
    signal +
    Array(probe.values) +
    Array(broadcast.values)
  }
  
  init() {
    for layerID in 1...4 {
      let y = 6 * Float(layerID)
      
      // Create 'signal'.
      do {
        let offset = SIMD3(0, y, 30.75)
        let pattern = PropagateUnit
          .signalPattern(layerID: layerID)
        let rod = PropagateUnit
          .createRodX(offset: offset, pattern: pattern)
        signal.append(rod)
      }
      
      // Create 'broadcast'.
      for positionX in 0..<layerID {
        var offset: SIMD3<Float>
        var pattern: KnobPattern
        
        if layerID == 4 && positionX == 3 {
          // Stack the final broadcast on the top layer, removing a large
          // block of unnecessary housing.
          let x = 7.5 * Float(positionX)
          offset = SIMD3(x + 11, y + 2.75, 0)
          pattern = { h, h2k, l in
            Concave {
              Convex {
                Origin { 46 * h }
                Plane { h }
              }
              Convex {
                Origin { 0.5 * h2k }
                Plane { -h2k }
              }
              Convex {
                Origin { 51 * h }
                Plane { -h }
              }
              Replace { .empty }
            }
          }
        } else {
          let x = 7.5 * Float(positionX)
          offset = SIMD3(x + 16, y - 2.75, 0)
          pattern = PropagateUnit
            .broadcastPattern()
        }
        let rod = PropagateUnit
          .createRodZ(offset: offset, pattern: pattern)
        
        let key = SIMD2(Int(positionX), Int(layerID))
        broadcast[key] = rod
      }
    }
    
    // Create 'probe'.
    for positionX in 0..<3 {
      let x = 7.5 * Float(positionX)
      let offset = SIMD3(x + 13.5, 0, 28)
      let pattern = PropagateUnit
        .probePattern(positionX: positionX)
      let rod = PropagateUnit
        .createRodY(offset: offset, pattern: pattern)
      
      let key = positionX
      probe[key] = rod
    }
  }
}

extension PropagateUnit {
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
      position = SIMD3(position.z, position.x, position.y)
      position += SIMD3(0.91, 0, 0.85)
      position += offset * latticeConstant
      copy.position = position
      return copy
    }
    return Rod(atoms: atoms)
  }
  
  private static func createRodZ(
    offset: SIMD3<Float>,
    pattern: KnobPattern
  ) -> Rod {
    let rodLatticeZ = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 54 * h + 2 * h2k + 2 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        pattern(h, h2k, l)
      }
    }
    
    let atoms = rodLatticeZ.atoms.map {
      var copy = $0
      var position = copy.position
      let latticeConstant = Constant(.square) { .elemental(.carbon) }
      position = SIMD3(position.z, position.y, position.x)
      position += SIMD3(0.91, 0.85, 0)
      position += offset * latticeConstant
      copy.position = position
      return copy
    }
    return Rod(atoms: atoms)
  }
}

extension PropagateUnit {
  private static func signalPattern(layerID: Int) -> KnobPattern {
    { h, h2k, l in
      let clockingShift: Float = 4
      
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
      
      // Create a groove for interaction with 'probe'.
      do {
        let startOffset: Float = 17 + clockingShift
        var endOffset: Float
        switch layerID {
        case 1: endOffset = 23 + clockingShift
        case 2: endOffset = 33 + clockingShift
        case 3: endOffset = 44 + clockingShift
        case 4: endOffset = 48 + clockingShift
        default: fatalError("Unexpected layer ID.")
        }
        
        Volume {
          Concave {
            Convex {
              Origin { startOffset * h }
              Plane { h }
            }
            Convex {
              Origin { 0.5 * l }
              Plane { -l }
            }
            Convex {
              Origin { endOffset * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
        createPhosphorusDopant(position: SIMD3(startOffset + 0.0, 0.75, -0.15))
        createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 1.75, 0.65))
        createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 0.25, 0.65))
        createPhosphorusDopant(position: SIMD3(endOffset + 0.0, 0.75, -0.15))
        createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 1.75, 0.65))
        createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 0.25, 0.65))
      }
      
      // Create a groove to avoid interaction with 'probe' on other layers.
      if layerID <= 2 {
        var startOffset: Float
        switch layerID {
        case 1: startOffset = 27.5 + clockingShift
        case 2: startOffset = 38.5 + clockingShift
        default: fatalError("Unexpected layer ID.")
        }
        let endOffset: Float = 47.5 + clockingShift
        
        Volume {
          Concave {
            Convex {
              Origin { startOffset * h }
              Plane { h }
            }
            Convex {
              Origin { 0.5 * l }
              Plane { -l }
            }
            Convex {
              Origin { endOffset * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
        createPhosphorusDopant(position: SIMD3(startOffset, 1.25, -0.15))
        createPhosphorusDopant(position: SIMD3(startOffset, 0.25, 0.65))
        createPhosphorusDopant(position: SIMD3(endOffset, 1.25, -0.15))
        createPhosphorusDopant(position: SIMD3(endOffset, 0.25, 0.65))
      }
      
      // Create a groove to directly transmit signals to 'broadcast'.
      if layerID == 4 {
        Volume {
          Concave {
            Convex {
              Origin { (46 + clockingShift) * h }
              Plane { h }
            }
            Convex {
              Origin { 1.5 * h2k }
              Plane { h2k }
            }
            Convex {
              Origin { (51 + clockingShift) * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
      
      func createPhosphorusDopant(position: SIMD3<Float>) {
        Volume {
          Concave {
            Concave {
              Origin { (position.x - 0.5) * h }
              Plane { h }
              Origin { 1 * h }
              Plane { -h }
            }
            Concave {
              Origin { (position.y - 0.25) * h2k }
              Plane { h2k }
              Origin { 0.5 * h2k }
              Plane { -h2k }
            }
            Concave {
              Origin { (position.z - 0.25) * l }
              Plane { l }
              Origin { 0.5 * l }
              Plane { -l }
            }
            Replace { .atom(.carbon) }
          }
        }
      }
    }
  }
  
  private static func probePattern(positionX: Int) -> KnobPattern {
    { h, h2k, l in
      let clockingShift: Float = 4
      
      // Create a groove to receive signals from 'signal'.
      do {
        Volume {
          Concave {
            Convex {
              var origin: Float
              switch positionX {
              case 0: origin = 11
              case 1: origin = 19
              case 2: origin = 28
              default: fatalError("Unrecognized position X.")
              }
              Origin { origin * h }
              Plane { h }
            }
            Convex {
              Origin { 1.5 * h2k }
              Plane { h2k }
            }
            Convex {
              var origin: Float
              switch positionX {
              case 0: origin = 16
              case 1: origin = 25
              case 2: origin = 33
              default: fatalError("Unrecognized position X.")
              }
              Origin { origin * h }
              Plane { -h }
            }
            Replace { .empty }
          }
        }
      }
      
      // Create a groove to transmit signals to 'broadcast' on layer 1.
      if positionX <= 0 {
        let startOffset: Float = 6.5 - clockingShift
        let endOffset: Float = 12.5 - clockingShift
        
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
        createPhosphorusDopant(position: SIMD3(startOffset, 1.25, 1.85))
        createPhosphorusDopant(position: SIMD3(startOffset, 0.25, 1.05))
        createPhosphorusDopant(position: SIMD3(endOffset, 1.25, 1.85))
        createPhosphorusDopant(position: SIMD3(endOffset, 0.25, 1.05))
      }
      
      // Create a groove to transmit signals to 'broadcast' on layer 2.
      if positionX <= 1 {
        let startOffset: Float = 15 - clockingShift
        let endOffset: Float = 21 - clockingShift
        
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
        createPhosphorusDopant(position: SIMD3(startOffset + 0.0, 0.75, 1.85))
        if positionX == 0 {
          createPhosphorusDopant(position: SIMD3(startOffset + 0.0, 1.75, 1.05))
        } else {
          createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 1.75, 1.05))
          createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 0.25, 1.05))
        }
        createPhosphorusDopant(position: SIMD3(endOffset + 0.0, 0.75, 1.85))
        createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 1.75, 1.05))
        if positionX != 0 {
          createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 0.25, 1.05))
        }
      }
      
      // Create a groove to transmit signals to 'broadcast' on layer 3.
      do {
        let startOffset: Float = 23.5 - clockingShift
        let endOffset: Float = 29.5 - clockingShift
        
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
        if positionX == 1 {
          createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 1.25, 1.05))
        } else {
          createPhosphorusDopant(position: SIMD3(startOffset, 1.25, 1.85))
          createPhosphorusDopant(position: SIMD3(startOffset, 0.25, 1.05))
        }
        createPhosphorusDopant(position: SIMD3(endOffset, 1.25, 1.85))
        createPhosphorusDopant(position: SIMD3(endOffset, 0.25, 1.05))
      }
      
      // Create a groove to transmit signals to 'broadcast' on layer 4.
      do {
        let startOffset: Float = 32 - clockingShift
        let endOffset: Float = 38 - clockingShift
        
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
        createPhosphorusDopant(position: SIMD3(startOffset + 0.0, 0.75, 1.85))
        if positionX == 2 {
          createPhosphorusDopant(position: SIMD3(startOffset + 0.0, 1.75, 1.05))
        } else {
          createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 1.75, 1.05))
          createPhosphorusDopant(position: SIMD3(startOffset + 0.5, 0.25, 1.05))
        }
        createPhosphorusDopant(position: SIMD3(endOffset + 0.0, 0.75, 1.85))
        createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 1.75, 1.05))
        createPhosphorusDopant(position: SIMD3(endOffset - 0.5, 0.25, 1.05))
      }
      
      func createPhosphorusDopant(position: SIMD3<Float>) {
        Volume {
          Concave {
            Concave {
              Origin { (position.x - 0.5) * h }
              Plane { h }
              Origin { 1 * h }
              Plane { -h }
            }
            Concave {
              Origin { (position.y - 0.25) * h2k }
              Plane { h2k }
              Origin { 0.5 * h2k }
              Plane { -h2k }
            }
            Concave {
              Origin { (position.z - 0.25) * l }
              Plane { l }
              Origin { 0.5 * l }
              Plane { -l }
            }
            Replace { .atom(.phosphorus) }
          }
        }
      }
    }
  }
  
  // We haven't reached the level of detail where individual broadcasts get
  // unique patterns.
  private static func broadcastPattern() -> KnobPattern {
    { h, h2k, l in
      // Create a groove to avoid interaction with 'signal'.
      Volume {
        Concave {
          Convex {
            Origin { 45 * h }
            Plane { h }
          }
          Convex {
            Origin { 1.5 * h2k }
            Plane { h2k }
          }
          Replace { .empty }
        }
      }
      
      // Create a groove to receive signals from 'probe'.
      Volume {
        Concave {
          Convex {
            Origin { 42 * h }
            Plane { h }
          }
          Convex {
            Origin { 0.5 * l }
            Plane { -l }
          }
          Convex {
            Origin { 48 * h }
            Plane { -h }
          }
          Replace { .empty }
        }
      }
      createPhosphorusDopant(position: SIMD3(42.0, 0.75, -0.15))
      createPhosphorusDopant(position: SIMD3(42.5, 1.75, 0.65))
      createPhosphorusDopant(position: SIMD3(42.5, 0.25, 0.65))
      createPhosphorusDopant(position: SIMD3(48.0, 0.75, -0.15))
      createPhosphorusDopant(position: SIMD3(47.5, 1.75, 0.65))
      createPhosphorusDopant(position: SIMD3(47.5, 0.25, 0.65))
      
      func createPhosphorusDopant(position: SIMD3<Float>) {
        Volume {
          Concave {
            Concave {
              Origin { (position.x - 0.5) * h }
              Plane { h }
              Origin { 1 * h }
              Plane { -h }
            }
            Concave {
              Origin { (position.y - 0.25) * h2k }
              Plane { h2k }
              Origin { 0.5 * h2k }
              Plane { -h2k }
            }
            Concave {
              Origin { (position.z - 0.25) * l }
              Plane { l }
              Origin { 0.5 * l }
              Plane { -l }
            }
            Replace { .atom(.phosphorus) }
          }
        }
      }
    }
  }
}
