//
//  BackBoard.swift
//  MolecularRendererApp
//
//  Created by Philip Turner on 2/24/24.
//

import Foundation
import HDL
import MM4
import Numerics
import OpenMM

struct BackBoard {
  var topology = Topology()
  
  init() {
    let backBoardLattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 400 * h + 100 * h2k + 6 * l }
      Material { .checkerboard(.silicon, .carbon) }
      
      Volume {
        Origin { 90 * h + 50 * h2k }
        
        func createHexagon() {
          var hexagonDirections: [SIMD3<Float>] = []
          hexagonDirections.append(k + 2 * h)
          hexagonDirections.append(h + 2 * k)
          hexagonDirections.append(k - h)
          hexagonDirections += hexagonDirections.map(-)
          
          for directionID in hexagonDirections.indices {
            let mainDirection = hexagonDirections[directionID]
            Convex {
              Origin { 43 * mainDirection }
              Plane { mainDirection }
            }
            Concave {
              let direction2 = hexagonDirections[(directionID + 1) % 6]
              let direction3 = hexagonDirections[(directionID + 5) % 6]
              Concave {
                Origin { 35 * mainDirection }
                Plane { -mainDirection }
                Convex {
                  Origin { 10 * direction2 }
                  Plane { -direction2 }
                }
                Convex {
                  Origin { 10 * direction3 }
                  Plane { -direction3 }
                }
              }
              Convex {
                Origin { 4 * direction2 }
                Plane { direction2 }
              }
              Convex {
                Origin { 4 * direction3 }
                Plane { direction3 }
              }
              Convex {
                Origin { 11.5 * mainDirection }
                Plane { mainDirection }
              }
            }
          }
        }
        
        Concave {
          Convex {
            createHexagon()
          }
          
          Convex {
            Plane { -h }
            Convex {
              Origin { 4 * h2k }
              Plane { h2k }
            }
            Convex {
              Origin { -4 * h2k }
              Plane { -h2k }
            }
          }
          
          Convex {
            Plane { -h2k }
            Convex {
              Origin { 86 * (k + h) }
              Plane { h2k }
              Plane { k + h + h2k / 3 }
            }
            
            Concave {
              Convex {
                Origin { 86 * (k + h) }
                Origin { -9 * (k + h + h2k / 3) }
                Plane { -(k + h + h2k / 3) }
              }
              Convex {
                Origin { 86 * (k + h) }
                Plane { -k - 2 * h }
                Origin { -30 * k }
                Plane { -k + h }
              }
              Convex {
                Origin { 171 * h }
                Plane { k - h }
              }
              Convex {
                Origin { 183 * h }
                Plane { -k - 2 * h }
              }
            }
          }
          Convex {
            Plane { h2k }
            Convex {
              Origin { 86 * (-k) }
              Plane { -h2k }
              Plane { -k - h2k / 3 }
            }
            
            Concave {
              Convex {
                Origin { 86 * (-k) }
                Origin { -9 * (-k - h2k / 3) }
                Plane { -(-k - h2k / 3) }
              }
              Convex {
                Origin { 86 * (-k) }
                Plane { k - h }
                Origin { -30 * (-k - h) }
                Plane { k + 2 * h }
              }
              Convex {
                Origin { 171 * h }
                Plane { -k - 2 * h }
              }
              Convex {
                Origin { 183 * h }
                Plane { k - h }
              }
            }
          }
        }
        
        Replace { .empty }
      }
    }
    
    var backBoardAtoms = backBoardLattice.atoms
    var accumulator: SIMD3<Double> = .zero
    var mass: Double = .zero
    for atom in backBoardAtoms {
      let atomMass = Float(atom.atomicNumber)
      accumulator += SIMD3<Double>(atomMass * atom.position)
      mass += Double(atomMass)
    }
    let centerOfMass = SIMD3<Float>(accumulator / mass)
    for i in backBoardAtoms.indices {
      backBoardAtoms[i].position -= centerOfMass
    }
    
    topology.insert(atoms: backBoardAtoms)
    var reconstruction = SurfaceReconstruction()
    reconstruction.material = .checkerboard(.silicon, .carbon)
    reconstruction.topology = topology
    reconstruction.removePathologicalAtoms()
    reconstruction.createBulkAtomBonds()
    reconstruction.createHydrogenSites()
    reconstruction.resolveCollisions()
    reconstruction.createHydrogenBonds()
    topology = reconstruction.topology
    topology.sort()
  }
}
