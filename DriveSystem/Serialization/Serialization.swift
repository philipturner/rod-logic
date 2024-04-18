//
//  Serialization.swift
//  MolecularRendererApp
//
//  Created by Philip Turner on 4/12/24.
//

import Foundation
import HDL
import MM4
import Numerics

struct Serialization {
  // Recover the bond topology after serializing the atom positions.
  //
  // The search algorithm uses a covalent bond scale of 1.4, instead of the
  // value of 1.5 from the HDL.
  static func reconstructBonds(
    topology: inout Topology,
    quaternaryAtomIDs: [UInt8]
  ) {
    // Only accept the 5 closest neighbors in the list.
    let matches = topology.match(
      topology.atoms,
      algorithm: .covalentBondLength(1.4),
      maximumNeighborCount: 5)
    
    var insertedBonds: [SIMD2<UInt32>] = []
    for atomID in topology.atoms.indices {
      let atom = topology.atoms[atomID]
      let neighborIDs = matches[atomID]
      if atom.atomicNumber == 1 {
        if neighborIDs.count != 2 {
          fatalError("Unexpected bond count: \(neighborIDs.count) != 2")
        }
      }
      if quaternaryAtomIDs.contains(atom.atomicNumber) {
        if neighborIDs.count < 5 {
          fatalError("Unexpected bond count: \(neighborIDs.count) < 5")
        }
      }
      
      for neighborID in neighborIDs where atomID < neighborID {
        let bond = SIMD2(UInt32(atomID), UInt32(neighborID))
        insertedBonds.append(bond)
      }
    }
    topology.insert(bonds: insertedBonds)
  }
  
  // Stores the raw data as a base-64 string.
  //
  // Compresses the information to within a factor of ~2x of the
  // information-theoretic limit.
  static func serialize(atoms: [Entity]) -> Data {
    // Allocate an array for the raw data.
    var rawData: [SIMD4<UInt16>] = []
    
    // Write to the header.
    var header: SIMD4<UInt16> = .zero
    header[0] = UInt16(atoms.count % 65536)
    header[1] = UInt16(atoms.count / 65536)
    rawData.append(header)
    
    for atom in atoms {
      // Quantize the position to ~1 picometer precision.
      var position = atom.position
      position *= 1024
      position += 32768
      position.round(.toNearestOrEven)
      
      // Check that every atom falls within the bounds: 0 ≤ r < 64 nm.
      if any(position .< 0) || any(position .>= 65536) {
        fatalError("Position was out of range.")
      }
      
      let chunk = SIMD4(
        SIMD3<UInt16>(position), UInt16(atom.atomicNumber))
      rawData.append(chunk)
    }
    
    // Return the data for further processing in the caller.
    return Data(bytes: rawData, count: rawData.count * 8)
  }
  
  // Retrieves the raw data from a base-64 string.
  static func deserialize(atoms: Data?) -> [Entity] {
    guard let data = atoms else {
      fatalError("Could not decode base64-encoded string.")
    }
    guard data.count % 8 == 0 else {
      fatalError("Data was not aligned to 8 bytes.")
    }
    
    var rawData = [SIMD4<UInt16>](repeating: .zero, count: data.count / 8)
    rawData.withUnsafeMutableBytes {
      let copiedCount = data.copyBytes(to: $0)
      guard copiedCount == data.count else {
        fatalError("Incorrect number of bytes was copied.")
      }
    }
    
    // Read the header.
    let header = rawData[0]
    let atomCount = UInt32(header[0]) + 65536 * UInt32(header[1])
    guard atomCount == rawData.count - 1 else {
      fatalError("Atom count was incorrect.")
    }
    
    var atoms: [Entity] = []
    for chunk in rawData[1...] {
      var position = SIMD3<Float>(
        unsafeBitCast(chunk, to: SIMD3<UInt16>.self))
      position -= 32768
      position /= 1024
      
      let atomicNumber = UInt8(chunk[3])
      guard let element = Element(rawValue: atomicNumber) else {
        fatalError("Unrecognized atomic number.")
      }
      let entity = Entity(position: position, type: .atom(element))
      atoms.append(entity)
    }
    return atoms
  }
  
  // Stores the raw data as a base-64 string.
  //
  // The typical compression ratio is ~2.5x. This is close to the
  // information-theoretic limit of ~3.5x.
  static func serialize(bonds: [SIMD2<UInt32>]) -> Data {
    var lastAtomID: UInt32 = .zero
    
    // Sort the bonds into compressable and non-compressable groups.
    var compressedBonds: [SIMD2<UInt8>] = []
    var decompressedBonds: [SIMD2<UInt32>] = []
    for bond in bonds {
      // We rely on there always being an 8-compressed atom nearby.
      guard bond[0] >= lastAtomID,
            bond[0] - lastAtomID < 256 else {
        fatalError("Bond could not be compressed.")
      }
      guard bond[1] > bond[0] else {
        fatalError("Bond was not sorted.")
      }
      
      let difference = bond[1] - bond[0]
      if difference < 256 {
        let startDelta = UInt8(bond[0] - lastAtomID)
        let lengthDelta = UInt8(bond[1] - bond[0])
        let compressedBond = SIMD2(startDelta, lengthDelta)
        compressedBonds.append(compressedBond)
        
        // Update the atom cursor.
        lastAtomID = bond[0]
      } else {
        // Do not update the atom cursor.
        decompressedBonds.append(bond)
      }
    }
    
    // Allocate an array for the raw data.
    var rawData: [SIMD2<UInt32>] = []
    
    // Write to the header.
    var header: SIMD2<UInt32> = .zero
    header[0] = UInt32(compressedBonds.count)
    header[1] = UInt32(decompressedBonds.count)
    rawData.append(header)
    
    // Pad the compressable bonds to a multiple of four.
    while compressedBonds.count % 4 != 0 {
      compressedBonds.append(SIMD2<UInt8>.zero)
    }
    
    // Write the compressable bonds.
    for groupID in 0..<compressedBonds.count / 4 {
      var vector: SIMD8<UInt8> = .zero
      for laneID in 0..<4 {
        let compressedBond = compressedBonds[groupID * 4 + laneID]
        vector[laneID * 2 + 0] = compressedBond[0]
        vector[laneID * 2 + 1] = compressedBond[1]
      }
      let castedVector = unsafeBitCast(vector, to: SIMD2<UInt32>.self)
      rawData.append(castedVector)
    }
    
    // Write the non-compressable bonds.
    rawData.append(contentsOf: decompressedBonds)
    
    // Return the data for further processing in the caller.
    return Data(bytes: rawData, count: rawData.count * 8)
  }
  
  // Retrieves the raw data from a base-64 string.
  static func deserialize(bonds: Data?) -> [SIMD2<UInt32>] {
    guard let data = bonds else {
      fatalError("Could not decode base64-encoded string.")
    }
    guard data.count % 8 == 0 else {
      fatalError("Data was not aligned to 8 bytes.")
    }
    
    var rawData = [SIMD2<UInt32>](repeating: .zero, count: data.count / 8)
    rawData.withUnsafeMutableBytes {
      let copiedCount = data.copyBytes(to: $0)
      guard copiedCount == data.count else {
        fatalError("Incorrect number of bytes was copied.")
      }
    }
    
    // Read the header.
    let header = rawData[0]
    let compressedBondCount = UInt32(header[0])
    let decompressedBondCount = UInt32(header[1])
    
    // Pad the compressable bonds to a multiple of four.
    var paddedCompressedBondCount = compressedBondCount
    while paddedCompressedBondCount % 4 != 0 {
      paddedCompressedBondCount += 1
    }
    
    // Read the compressable bonds.
    var compressedBonds: [SIMD2<UInt8>] = []
    for groupID in 0..<paddedCompressedBondCount / 4 {
      let castedVector = rawData[Int(1 + groupID)]
      let vector = unsafeBitCast(castedVector, to: SIMD8<UInt8>.self)
      for laneID in 0..<4 {
        var compressedBond: SIMD2<UInt8> = .zero
        compressedBond[0] = vector[laneID * 2 + 0]
        compressedBond[1] = vector[laneID * 2 + 1]
        compressedBonds.append(compressedBond)
      }
    }
    
    // Restore the padded count to the original count.
    while compressedBonds.count > compressedBondCount {
      compressedBonds.removeLast()
    }
    
    // Read the non-compressable bonds.
    let decompressedRange = Int(1 + paddedCompressedBondCount / 4)...
    let decompressedBonds = Array(rawData[decompressedRange])
    guard decompressedBonds.count == decompressedBondCount else {
      fatalError("Could not decode decompressed bonds.")
    }
    
    // Merge the bonds into a common array.
    var lastAtomID: UInt32 = .zero
    var bonds: [SIMD2<UInt32>] = []
    for compressedBond in compressedBonds {
      let startDelta = UInt32(compressedBond[0])
      let lengthDelta = UInt32(compressedBond[1])
      let bond = SIMD2(
        lastAtomID + startDelta,
        lastAtomID + startDelta + lengthDelta)
      bonds.append(bond)
      
      // Update the atom cursor.
      lastAtomID += startDelta
    }
    bonds.append(contentsOf: decompressedBonds)
    
    // Sort the array of bonds.
    bonds.sort(by: { lhs, rhs in
      if lhs[0] > rhs[0] {
        return false
      } else if lhs[0] < rhs[0] {
        return true
      }
      if lhs[1] > rhs[1] {
        return false
      } else if lhs[1] < rhs[1] {
        return true
      }
      return true
    })
    return bonds
  }
}
