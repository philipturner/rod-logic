// Save as GitHub Gist instead of polluting the molecular-renderer source tree.

import Foundation
import HDL
import MM4
import Numerics

func createGeometry() -> [[Entity]] {
  var output: [Entity] = []
  let useNitrogen = true
  let housing = createHousing(useNitrogen: useNitrogen)
  let rod = createRod(useNitrogen: useNitrogen).map {
    var copy = $0
//    copy.position.y += 0.4
    copy.position.z -= 1
    return copy
  }
  
  var fluorines = housing.filter { $0.atomicNumber == 9 }
  var silicons = housing.filter { $0.atomicNumber == 14 }
  let carbons = housing.filter { $0.atomicNumber == 6 }
  let hydrogens = housing.filter { $0.atomicNumber == 1 }
  print("GFN2-xTB atoms (compute cost):", 2 * fluorines.count + silicons.count + rod.count)
  print("GFN2-xTB atoms:", fluorines.count + silicons.count + rod.count)
  print("GFN-FF atoms:", housing.count + rod.count)
  
  for i in fluorines.indices {
    fluorines[i].atomicNumber = 1
  }
  for i in silicons.indices {
    silicons[i].atomicNumber = 6
  }
  
  var topology1 = Topology()
  topology1.insert(atoms: carbons + hydrogens)
  topology1.sort()
  var topology2 = Topology()
  topology2.insert(atoms: fluorines + silicons)
  topology2.sort()
  
  output = topology1.atoms + topology2.atoms + rod
  let expensiveAtomStart = topology1.atoms.count
  
  
  setenv("OMP_STACKSIZE", "2G", 1)
  setenv("OMP_NUM_THREADS", "7", 1)
  let path = "/Users/philipturner/Documents/OpenMM/xtb/cpu0"
  let process = XTBProcess(path: path)
  process.engine = "inertial"
  process.writeFile(name: "xtb.inp", process.encodeSettings())
//  process.writeFile(name: "coord", try! process.encodeAtoms(output))
//  process.run(arguments: [
//    "coord", "--gfnff", "--opt", "crude", "--input", "xtb.inp",
//  ])
//  output = try! process.decodeAtoms(process.readFile(name: "xtbopt.coord"))
  
  let innerRegion = "\(expensiveAtomStart + 1)-\(output.count)"
  print(expensiveAtomStart)
  print(innerRegion)
//  process.writeFile(name: "coord", try! process.encodeAtoms(output))
//  process.run(arguments: ["coord", "--oniom", "gfn2:gfnff", innerRegion, "--opt", "crude", "--input", "xtb.inp"])
//  output = try! process.decodeAtoms(process.readFile(name: "xtbopt.coord"))
//  
  // Next steps:
  // - test diamond system in GFN-FF
  // - test nitrogen system in GFN-FF
  // - test nitrogen system in ONIOM
  
  for i in output.indices {
    output[i].position += SIMD3(-1.35, -1.25, -2)
  }
  return [output, Array(output[expensiveAtomStart...])]
}

// MARK: - Geometry Compilation

func createHousing(useNitrogen: Bool) -> [Entity] {
  let housingLattice = Lattice<Cubic> { h, k, l in
    Bounds { 8 * h + 5 * k + 5 * l }
    Material { .elemental(.carbon) }
    
    Volume {
      Concave {
        let planePosition = Float(6) - (useNitrogen ? 0.5 : 0)
        Convex {
          Origin { 1.5 * k }
          Plane { k }
        }
        Convex {
          Origin { 2 * h }
          Plane { h }
        }
        Convex {
          Origin { planePosition * h }
          Plane { -h }
        }
        Convex {
          Origin { 2 * h + 2 * k }
          Plane { h + k }
        }
        Convex {
          Origin { planePosition * h + 2 * k }
          Plane { -h + k }
        }
      }
      
      let planePosition = Float(7.5) - (useNitrogen ? 0.5 : 0)
      Convex {
        Origin { 0.5 * h }
        Plane { -h }
      }
      Convex {
        Origin { planePosition * h }
        Plane { h }
      }
      Convex {
        Origin { 0.5 * l }
        Plane { -l }
      }
      Convex {
        Origin { 4.5 * l }
        Plane { l }
      }
      
      Replace { .empty }
      
      Concave {
        let planePosition = Float(6.25) - (useNitrogen ? 0.5 : 0)
        Convex {
          Origin { 1.25 * k }
          Plane { k }
        }
        Convex {
          Origin { 1.75 * h }
          Plane { h }
        }
        Convex {
          Origin { planePosition * h }
          Plane { -h }
        }
        Convex {
          Origin { 1.75 * h + 1.75 * k }
          Plane { h + k }
        }
        Convex {
          Origin { planePosition * h + 1.75 * k }
          Plane { -h + k }
        }
        Convex {
          Origin { 0.5 * l }
          Plane { l }
        }
        Convex {
          Origin { 4.5 * l }
          Plane { -l }
        }
      }
      
      Replace { .atom(.silicon) }
    }
  }
  
  var topology = Topology()
  topology.insert(atoms: housingLattice.atoms)
  
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
  
  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = topology.atoms.map(\.atomicNumber)
  for i in paramsDesc.atomicNumbers!.indices {
    if paramsDesc.atomicNumbers![i] == 14 {
      paramsDesc.atomicNumbers![i] = 6
    }
  }
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
  
  let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
  for i in topology.atoms.indices {
    let neighbors = atomsToAtomsMap[i]
    if topology.atoms[i].atomicNumber == 1 {
      let neighborID = Int(neighbors.first!)
      if topology.atoms[neighborID].atomicNumber == 14 {
        topology.atoms[i].atomicNumber = 9
      }
    }
  }
  
  return topology.atoms
}

func createRod(useNitrogen: Bool) -> [Entity] {
  let rodLattice = Lattice<Hexagonal> { h, k, l in
    let h2k = h + 2 * k
    Bounds { Float(8) * h + 2 * h2k + 4 * l }
    Material { .elemental(.carbon) }
    
    Volume {
      Convex {
        Origin { 1.9 * l }
        Plane { l }
      }
      Concave {
        Convex {
          Origin { Float(1) * h }
          Plane { h }
        }
        Convex {
          Origin { 1.5 * h2k }
          Plane { h2k }
        }
        Convex {
          Origin { Float(10) * h }
          Plane { -h }
        }
      }
      Replace { .empty }
    }
  }
  
  var topology = Topology()
  topology.insert(atoms: rodLattice.atoms)
  
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
  
  if useNitrogen {
    let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
    var removedAtoms: [UInt32] = []
    for i in topology.atoms.indices {
      let position = topology.atoms[i].position
      if position.z < 0.05 || position.z > 0.65 {
        if topology.atoms[i].atomicNumber == 6 {
          topology.atoms[i].atomicNumber = 7
        } else {
          
        }
      }
    }
    for i in topology.atoms.indices {
      let position = topology.atoms[i].position
      if position.z < 0.00 || position.z > 0.70 {
        if topology.atoms[i].atomicNumber == 6 {
          
        } else {
          let neighborID = Int(atomsToAtomsMap[i].first!)
          if topology.atoms[neighborID].atomicNumber == 7 {
            removedAtoms.append(UInt32(i))
          }
        }
      }
    }
    topology.remove(atoms: removedAtoms)
  }
  
  let latticeConstant = Constant(.square) { .elemental(.carbon) }
  var translation: SIMD3<Float> = .zero
  translation.x += 0.025 + 3 * latticeConstant
  translation.y += -0.050 + 2.5 * latticeConstant
  translation.z += 1.0 * latticeConstant
  if useNitrogen {
    translation.x -= 0.25 * latticeConstant
  }
  
  return topology.atoms.map {
    var copy = $0
    copy.position = SIMD3(copy.position.z, copy.position.y, copy.position.x)
    copy.position += translation
    return copy
  }
}

// MARK: - xTB Outputs

/*
GFN2-xTB atoms (compute cost): 230
GFN2-xTB atoms: 196
GFN-FF atoms: 832
636
637-832
*/

// Rod inside housing:

/*
 ........................................................................
 .............................. CYCLE   96 ..............................
 ........................................................................

   ------------------------------------------------------------------------

       Singlepoint calculation of whole system with low-level method

   ------------------------------------------------------------------------


  E+G (total)                   0 d,  0 h,  0 min,  0.101 sec
  distance/D3 list               ...        0 min,  0.001 sec (  0.851%)
  non bonded repulsion           ...        0 min,  0.001 sec (  1.169%)
  dCN                            ...        0 min,  0.003 sec (  3.084%)
  EEQ energy and q               ...        0 min,  0.039 sec ( 38.483%)
  D3                             ...        0 min,  0.051 sec ( 50.550%)
  EEQ gradient                   ...        0 min,  0.001 sec (  1.088%)
  bonds                          ...        0 min,  0.004 sec (  3.800%)
  bend and torsion               ...        0 min,  0.000 sec (  0.334%)
  bonded ATM                     ...        0 min,  0.000 sec (  0.288%)
  HB/XB (incl list setup)        ...        0 min,  0.000 sec (  0.328%)

          :: total energy            -162.041974221135 Eh    ::
          :: gradient norm              0.993214593789 Eh/a0 ::

   ------------------------------------------------------------------------

       Singlepoint calculation of inner region with low-level method

   ------------------------------------------------------------------------


  E+G (total)                   0 d,  0 h,  0 min,  0.011 sec
  distance/D3 list               ...        0 min,  0.000 sec (  0.806%)
  non bonded repulsion           ...        0 min,  0.000 sec (  2.314%)
  dCN                            ...        0 min,  0.001 sec (  7.480%)
  EEQ energy and q               ...        0 min,  0.002 sec ( 21.685%)
  D3                             ...        0 min,  0.006 sec ( 55.410%)
  EEQ gradient                   ...        0 min,  0.000 sec (  2.362%)
  bonds                          ...        0 min,  0.000 sec (  4.131%)
  bend and torsion               ...        0 min,  0.000 sec (  1.603%)
  bonded ATM                     ...        0 min,  0.000 sec (  0.870%)
  HB/XB (incl list setup)        ...        0 min,  0.000 sec (  3.148%)

          :: total energy             -37.608739606057 Eh    ::
          :: gradient norm              0.929836258917 Eh/a0 ::

   ------------------------------------------------------------------------

       Singlepoint calculation of inner region with high-level method

   ------------------------------------------------------------------------


  iter      E             dE          RMSdq      gap      omega  full diag
    1   -344.5949031 -0.344595E+03  0.177E-03    5.36       0.0  T
    2   -344.5949032 -0.903743E-07  0.102E-03    5.36       5.9  T
    3   -344.5949032  0.125027E-08  0.343E-04    5.36      17.5  T
    4   -344.5949032 -0.822894E-08  0.101E-04    5.36      59.3  T
      SCC iter.                  ...        0 min,  0.623 sec
      gradient                   ...        0 min,  0.271 sec
  * total energy  :  -464.4323912 Eh     change    :  -0.3201771E-03 Eh
    gradient norm :     0.0085001 Eh/a0  converged?    E=T G=T D=F
    time step     :     0.2500000 fs     power     :   0.1655647E-07
    step length   :     0.0613984 a0     speed     :   0.5940625E-02

    *** GEOMETRY OPTIMIZATION CONVERGED AFTER 96 CYCLES ***

 ------------------------------------------------------------------------
  total energy gain   :        -1.9961686 Eh    -1252.6147 kcal/mol
  total RMSD          :         0.3625855 a0        0.1919 Å
  total power (kW/mol):       -54.5931238 (step)  -35.7501 (real)
 ------------------------------------------------------------------------

  FIRE (total)                  0 d,  0 h,  2 min, 26.600 sec
  setup                          ...        0 min,  0.003 sec (  0.002%)
  velocities                     ...        0 min,  0.000 sec (  0.000%)
  propagation                    ...        0 min,  0.001 sec (  0.000%)
  singlepoint calculation        ...        2 min, 26.481 sec ( 99.919%)
  log and printout               ...        0 min,  0.112 sec (  0.076%)
  preconditioner                 ...        0 min,  0.000 sec (  0.000%)
  model hessian                  ...        0 min,  0.000 sec (  0.000%)
  hessian update                 ...        0 min,  0.000 sec (  0.000%)
 
 -------------------------------------------------
|                Final Singlepoint                |
 -------------------------------------------------

------------------------------------------------------------------------

Singlepoint calculation of whole system with low-level method

------------------------------------------------------------------------


E+G (total)                   0 d,  0 h,  0 min,  0.102 sec
distance/D3 list               ...        0 min,  0.001 sec (  0.921%)
non bonded repulsion           ...        0 min,  0.001 sec (  1.356%)
dCN                            ...        0 min,  0.003 sec (  3.296%)
EEQ energy and q               ...        0 min,  0.039 sec ( 38.263%)
D3                             ...        0 min,  0.052 sec ( 50.681%)
EEQ gradient                   ...        0 min,  0.001 sec (  1.209%)
bonds                          ...        0 min,  0.003 sec (  3.204%)
bend and torsion               ...        0 min,  0.000 sec (  0.355%)
bonded ATM                     ...        0 min,  0.000 sec (  0.280%)
HB/XB (incl list setup)        ...        0 min,  0.000 sec (  0.410%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy            -162.041974221135 Eh    ::
:: gradient norm              0.993214593789 Eh/a0 ::
::.................................................::
:: bond energy             -179.981748011284 Eh    ::
:: angle energy               1.895730049796 Eh    ::
:: torsion energy             2.661722218333 Eh    ::
:: repulsion energy          17.072365841911 Eh    ::
:: electrostat energy        -0.290193471426 Eh    ::
:: dispersion energy         -2.596937358865 Eh    ::
:: HB energy                 -0.007146519404 Eh    ::
:: XB energy                  0.000000000000 Eh    ::
:: bonded atm energy         -0.795766970196 Eh    ::
:: external energy            0.000000000000 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge              -0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::


------------------------------------------------------------------------

Singlepoint calculation of inner region with low-level method

------------------------------------------------------------------------


E+G (total)                   0 d,  0 h,  0 min,  0.011 sec
distance/D3 list               ...        0 min,  0.000 sec (  0.448%)
non bonded repulsion           ...        0 min,  0.000 sec (  2.464%)
dCN                            ...        0 min,  0.001 sec (  7.633%)
EEQ energy and q               ...        0 min,  0.003 sec ( 27.271%)
D3                             ...        0 min,  0.006 sec ( 52.133%)
EEQ gradient                   ...        0 min,  0.000 sec (  1.694%)
bonds                          ...        0 min,  0.000 sec (  3.410%)
bend and torsion               ...        0 min,  0.000 sec (  1.096%)
bonded ATM                     ...        0 min,  0.000 sec (  0.802%)
HB/XB (incl list setup)        ...        0 min,  0.000 sec (  2.907%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy             -37.608739606057 Eh    ::
:: gradient norm              0.929836258917 Eh/a0 ::
::.................................................::
:: bond energy              -41.731555780983 Eh    ::
:: angle energy               0.247759163675 Eh    ::
:: torsion energy             0.356771572378 Eh    ::
:: repulsion energy           4.384459345920 Eh    ::
:: electrostat energy        -0.221640553129 Eh    ::
:: dispersion energy         -0.491698356642 Eh    ::
:: HB energy                 -0.008044955387 Eh    ::
:: XB energy                  0.000000000000 Eh    ::
:: bonded atm energy         -0.144790041891 Eh    ::
:: external energy            0.000000000000 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge               0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::


------------------------------------------------------------------------

Singlepoint calculation of inner region with high-level method

------------------------------------------------------------------------


...................................................
:                      SETUP                      :
:.................................................:
:  # basis functions                 620          :
:  # atomic orbitals                 620          :
:  # shells                          392          :
:  # electrons                       638          :
:  max. iterations                   250          :
:  Hamiltonian                  GFN2-xTB          :
:  restarted?                      false          :
:  GBSA solvation                  false          :
:  PC potential                    false          :
:  electronic temp.          300.0000000     K    :
:  accuracy                    1.0000000          :
:  -> integral cutoff          0.2500000E+02      :
:  -> integral neglect         0.1000000E-07      :
:  -> SCF convergence          0.1000000E-05 Eh   :
:  -> wf. convergence          0.1000000E-03 e    :
:  Broyden damping             0.4000000          :
...................................................

iter      E             dE          RMSdq      gap      omega  full diag
1   -344.5949032 -0.344595E+03  0.326E-05    5.36       0.0  T
2   -344.5949032 -0.443379E-10  0.177E-05    5.36     338.1  T
3   -344.5949032  0.568434E-13  0.906E-06    5.36     661.9  T

*** convergence criteria satisfied after 3 iterations ***

#    Occupation            Energy/Eh            Energy/eV
-------------------------------------------------------------
1        2.0000           -0.6953960             -18.9227
...           ...                  ...                  ...
313        2.0000           -0.3273959              -8.9089
314        2.0000           -0.3254116              -8.8549
315        2.0000           -0.3235603              -8.8045
316        2.0000           -0.3173292              -8.6350
317        2.0000           -0.3161453              -8.6028
318        2.0000           -0.2901293              -7.8948
319        2.0000           -0.2891708              -7.8687 (HOMO)
320                         -0.0922219              -2.5095 (LUMO)
321                         -0.0883017              -2.4028
322                         -0.0792509              -2.1565
323                         -0.0761828              -2.0730
324                         -0.0576501              -1.5687
...                                ...                  ...
620                          1.1612478              31.5992
-------------------------------------------------------------
        HL-Gap            0.1969489 Eh            5.3593 eV
   Fermi-level           -0.1906964 Eh           -5.1891 eV

SCC (total)                   0 d,  0 h,  0 min,  0.789 sec
SCC setup                      ...        0 min,  0.001 sec (  0.144%)
Dispersion                     ...        0 min,  0.005 sec (  0.610%)
classical contributions        ...        0 min,  0.001 sec (  0.124%)
integral evaluation            ...        0 min,  0.020 sec (  2.546%)
iterations                     ...        0 min,  0.487 sec ( 61.734%)
molecular gradient             ...        0 min,  0.271 sec ( 34.290%)
printout                       ...        0 min,  0.004 sec (  0.549%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy            -339.999156565982 Eh    ::
:: gradient norm              1.002894360532 Eh/a0 ::
:: HOMO-LUMO gap              5.359252848273 eV    ::
::.................................................::
:: SCC energy              -344.594903246844 Eh    ::
:: -> isotropic ES            0.300654364905 Eh    ::
:: -> anisotropic ES          0.094065235615 Eh    ::
:: -> anisotropic XC          0.229789199838 Eh    ::
:: -> dispersion             -0.475444840164 Eh    ::
:: repulsion energy           4.530995690137 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge              -0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
 
 -------------------------------------------------
| TOTAL ENERGY             -464.432391181059 Eh   |
| GRADIENT NORM               0.008500429453 Eh/α |
 -------------------------------------------------

------------------------------------------------------------------------
* finished run on 2024/02/12 at 21:57:04.934
------------------------------------------------------------------------
total:
* wall-time:     0 d,  0 h,  2 min, 36.240 sec
*  cpu-time:     0 d,  0 h,  8 min, 19.884 sec
* ratio c/w:     3.199 speedup
SCF:
* wall-time:     0 d,  0 h,  0 min,  3.881 sec
*  cpu-time:     0 d,  0 h,  0 min, 14.403 sec
* ratio c/w:     3.711 speedup
ANC optimizer:
* wall-time:     0 d,  0 h,  2 min, 27.647 sec
*  cpu-time:     0 d,  0 h,  7 min, 42.342 sec
* ratio c/w:     3.131 speedup

normal termination of xtb
Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_UNDERFLOW_FLAG
atoms: 832
frames: 2
setup time: 156508.7 ms
 */

// Rod outside housing:

/*
 ........................................................................
 .............................. CYCLE   96 ..............................
 ........................................................................

   ------------------------------------------------------------------------

       Singlepoint calculation of whole system with low-level method

   ------------------------------------------------------------------------


  E+G (total)                   0 d,  0 h,  0 min,  0.104 sec
  distance/D3 list               ...        0 min,  0.001 sec (  0.861%)
  non bonded repulsion           ...        0 min,  0.001 sec (  0.825%)
  dCN                            ...        0 min,  0.003 sec (  2.763%)
  EEQ energy and q               ...        0 min,  0.047 sec ( 44.877%)
  D3                             ...        0 min,  0.047 sec ( 45.536%)
  EEQ gradient                   ...        0 min,  0.001 sec (  1.026%)
  bonds                          ...        0 min,  0.003 sec (  3.247%)
  bend and torsion               ...        0 min,  0.000 sec (  0.350%)
  bonded ATM                     ...        0 min,  0.000 sec (  0.256%)
  HB/XB (incl list setup)        ...        0 min,  0.000 sec (  0.241%)

          :: total energy            -161.915062344301 Eh    ::
          :: gradient norm              0.998496036159 Eh/a0 ::

   ------------------------------------------------------------------------

       Singlepoint calculation of inner region with low-level method

   ------------------------------------------------------------------------


  E+G (total)                   0 d,  0 h,  0 min,  0.010 sec
  distance/D3 list               ...        0 min,  0.000 sec (  0.988%)
  non bonded repulsion           ...        0 min,  0.000 sec (  1.840%)
  dCN                            ...        0 min,  0.001 sec (  6.228%)
  EEQ energy and q               ...        0 min,  0.003 sec ( 29.947%)
  D3                             ...        0 min,  0.005 sec ( 50.209%)
  EEQ gradient                   ...        0 min,  0.000 sec (  2.090%)
  bonds                          ...        0 min,  0.000 sec (  4.078%)
  bend and torsion               ...        0 min,  0.000 sec (  1.337%)
  bonded ATM                     ...        0 min,  0.000 sec (  0.824%)
  HB/XB (incl list setup)        ...        0 min,  0.000 sec (  2.265%)

          :: total energy             -37.543272576790 Eh    ::
          :: gradient norm              0.930199287435 Eh/a0 ::

   ------------------------------------------------------------------------

       Singlepoint calculation of inner region with high-level method

   ------------------------------------------------------------------------


  iter      E             dE          RMSdq      gap      omega  full diag
    1   -344.4954427 -0.344495E+03  0.168E-03    6.42       0.0  T
    2   -344.4954427 -0.693058E-07  0.972E-04    6.42       6.2  T
    3   -344.4954427  0.244972E-08  0.312E-04    6.42      19.2  T
      SCC iter.                  ...        0 min,  0.488 sec
      gradient                   ...        0 min,  0.172 sec
  * total energy  :  -464.3618554 Eh     change    :  -0.1905016E-03 Eh
    gradient norm :     0.0076317 Eh/a0  converged?    E=T G=T D=F
    time step     :     0.2500000 fs     power     :   0.9728961E-08
    step length   :     0.0523371 a0     speed     :   0.5063892E-02

    *** GEOMETRY OPTIMIZATION CONVERGED AFTER 96 CYCLES ***

 ------------------------------------------------------------------------
  total energy gain   :        -2.0097572 Eh    -1261.1417 kcal/mol
  total RMSD          :         0.3590094 a0        0.1900 Å
  total power (kW/mol):       -54.9647589 (step)  -45.0232 (real)
 ------------------------------------------------------------------------

  FIRE (total)                  0 d,  0 h,  1 min, 57.198 sec
  setup                          ...        0 min,  0.003 sec (  0.003%)
  velocities                     ...        0 min,  0.000 sec (  0.000%)
  propagation                    ...        0 min,  0.001 sec (  0.000%)
  singlepoint calculation        ...        1 min, 57.081 sec ( 99.900%)
  log and printout               ...        0 min,  0.111 sec (  0.094%)
  preconditioner                 ...        0 min,  0.000 sec (  0.000%)
  model hessian                  ...        0 min,  0.000 sec (  0.000%)
  hessian update                 ...        0 min,  0.000 sec (  0.000%)
 
 -------------------------------------------------
|                Final Singlepoint                |
 -------------------------------------------------

------------------------------------------------------------------------

Singlepoint calculation of whole system with low-level method

------------------------------------------------------------------------


E+G (total)                   0 d,  0 h,  0 min,  0.107 sec
distance/D3 list               ...        0 min,  0.001 sec (  0.959%)
non bonded repulsion           ...        0 min,  0.001 sec (  0.837%)
dCN                            ...        0 min,  0.003 sec (  2.874%)
EEQ energy and q               ...        0 min,  0.051 sec ( 47.543%)
D3                             ...        0 min,  0.045 sec ( 42.539%)
EEQ gradient                   ...        0 min,  0.001 sec (  0.924%)
bonds                          ...        0 min,  0.004 sec (  3.405%)
bend and torsion               ...        0 min,  0.000 sec (  0.347%)
bonded ATM                     ...        0 min,  0.000 sec (  0.290%)
HB/XB (incl list setup)        ...        0 min,  0.000 sec (  0.262%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy            -161.915062344301 Eh    ::
:: gradient norm              0.998496036159 Eh/a0 ::
::.................................................::
:: bond energy             -179.930830949784 Eh    ::
:: angle energy               1.900057987862 Eh    ::
:: torsion energy             2.650565570033 Eh    ::
:: repulsion energy          16.990336532689 Eh    ::
:: electrostat energy        -0.279338370473 Eh    ::
:: dispersion energy         -2.452958870084 Eh    ::
:: HB energy                 -0.000017665649 Eh    ::
:: XB energy                  0.000000000000 Eh    ::
:: bonded atm energy         -0.792876578896 Eh    ::
:: external energy            0.000000000000 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge               0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::


------------------------------------------------------------------------

Singlepoint calculation of inner region with low-level method

------------------------------------------------------------------------


E+G (total)                   0 d,  0 h,  0 min,  0.009 sec
distance/D3 list               ...        0 min,  0.000 sec (  0.907%)
non bonded repulsion           ...        0 min,  0.000 sec (  2.151%)
dCN                            ...        0 min,  0.001 sec (  7.196%)
EEQ energy and q               ...        0 min,  0.002 sec ( 25.980%)
D3                             ...        0 min,  0.005 sec ( 52.720%)
EEQ gradient                   ...        0 min,  0.000 sec (  1.816%)
bonds                          ...        0 min,  0.000 sec (  4.654%)
bend and torsion               ...        0 min,  0.000 sec (  1.332%)
bonded ATM                     ...        0 min,  0.000 sec (  0.919%)
HB/XB (incl list setup)        ...        0 min,  0.000 sec (  2.192%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy             -37.543272576790 Eh    ::
:: gradient norm              0.930199287435 Eh/a0 ::
::.................................................::
:: bond energy              -41.694846588118 Eh    ::
:: angle energy               0.249451270040 Eh    ::
:: torsion energy             0.349567507844 Eh    ::
:: repulsion energy           4.297855225119 Eh    ::
:: electrostat energy        -0.211516048525 Eh    ::
:: dispersion energy         -0.391710460762 Eh    ::
:: HB energy                 -0.000021989527 Eh    ::
:: XB energy                  0.000000000000 Eh    ::
:: bonded atm energy         -0.142051492860 Eh    ::
:: external energy            0.000000000000 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge               0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::


------------------------------------------------------------------------

Singlepoint calculation of inner region with high-level method

------------------------------------------------------------------------


...................................................
:                      SETUP                      :
:.................................................:
:  # basis functions                 620          :
:  # atomic orbitals                 620          :
:  # shells                          392          :
:  # electrons                       638          :
:  max. iterations                   250          :
:  Hamiltonian                  GFN2-xTB          :
:  restarted?                      false          :
:  GBSA solvation                  false          :
:  PC potential                    false          :
:  electronic temp.          300.0000000     K    :
:  accuracy                    1.0000000          :
:  -> integral cutoff          0.2500000E+02      :
:  -> integral neglect         0.1000000E-07      :
:  -> SCF convergence          0.1000000E-05 Eh   :
:  -> wf. convergence          0.1000000E-03 e    :
:  Broyden damping             0.4000000          :
...................................................

iter      E             dE          RMSdq      gap      omega  full diag
1   -344.4954427 -0.344495E+03  0.970E-05    6.42       0.0  T
2   -344.4954428 -0.100431E-08  0.474E-05    6.42     126.6  T
3   -344.4954428 -0.305533E-09  0.271E-05    6.42     221.2  T

*** convergence criteria satisfied after 3 iterations ***

#    Occupation            Energy/Eh            Energy/eV
-------------------------------------------------------------
1        2.0000           -0.6905564             -18.7910
...           ...                  ...                  ...
313        2.0000           -0.3466099              -9.4317
314        2.0000           -0.3444760              -9.3737
315        2.0000           -0.3412610              -9.2862
316        2.0000           -0.3269105              -8.8957
317        2.0000           -0.3256080              -8.8602
318        2.0000           -0.3216151              -8.7516
319        2.0000           -0.3195059              -8.6942 (HOMO)
320                         -0.0836559              -2.2764 (LUMO)
321                         -0.0807933              -2.1985
322                         -0.0707027              -1.9239
323                         -0.0688139              -1.8725
324                         -0.0668428              -1.8189
...                                ...                  ...
620                          1.1263714              30.6501
-------------------------------------------------------------
        HL-Gap            0.2358500 Eh            6.4178 eV
   Fermi-level           -0.2015809 Eh           -5.4853 eV

SCC (total)                   0 d,  0 h,  0 min,  0.654 sec
SCC setup                      ...        0 min,  0.001 sec (  0.157%)
Dispersion                     ...        0 min,  0.005 sec (  0.751%)
classical contributions        ...        0 min,  0.001 sec (  0.127%)
integral evaluation            ...        0 min,  0.019 sec (  2.882%)
iterations                     ...        0 min,  0.455 sec ( 69.531%)
molecular gradient             ...        0 min,  0.170 sec ( 26.000%)
printout                       ...        0 min,  0.004 sec (  0.548%)

:::::::::::::::::::::::::::::::::::::::::::::::::::::
::                     SUMMARY                     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::
:: total energy            -339.990065593082 Eh    ::
:: gradient norm              0.997945831572 Eh/a0 ::
:: HOMO-LUMO gap              6.417806487753 eV    ::
::.................................................::
:: SCC energy              -344.495442750625 Eh    ::
:: -> isotropic ES            0.294582057054 Eh    ::
:: -> anisotropic ES          0.126549200141 Eh    ::
:: -> anisotropic XC          0.232778198434 Eh    ::
:: -> dispersion             -0.401157794166 Eh    ::
:: repulsion energy           4.473664034654 Eh    ::
:: add. restraining           0.000000000000 Eh    ::
:: total charge              -0.000000000000 e     ::
:::::::::::::::::::::::::::::::::::::::::::::::::::::

 -------------------------------------------------
| TOTAL ENERGY             -464.361855360594 Eh   |
| GRADIENT NORM               0.007629734036 Eh/α |
 -------------------------------------------------

------------------------------------------------------------------------
* finished run on 2024/02/12 at 22:06:59.149
------------------------------------------------------------------------
total:
* wall-time:     0 d,  0 h,  2 min,  1.819 sec
*  cpu-time:     0 d,  0 h,  5 min, 28.916 sec
* ratio c/w:     2.700 speedup
SCF:
* wall-time:     0 d,  0 h,  0 min,  2.488 sec
*  cpu-time:     0 d,  0 h,  0 min,  6.225 sec
* ratio c/w:     2.502 speedup
ANC optimizer:
* wall-time:     0 d,  0 h,  1 min, 58.113 sec
*  cpu-time:     0 d,  0 h,  5 min, 20.353 sec
* ratio c/w:     2.712 speedup

normal termination of xtb
Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG
atoms: 832
frames: 2
setup time: 122047.4 ms
 */
