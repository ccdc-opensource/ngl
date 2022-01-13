/**
 * @file CCDCMol2 Parser
 * @author CCDC
 * @private
 */

import { Debug, Log, ParserRegistry } from '../globals'
import {
  assignResidueTypeBonds,
  calculateChainnames, calculateSecondaryStructure,
  calculateBondsBetween, calculateBondsWithin
} from '../structure/structure-utils'
import StructureParser from './structure-parser'
import { Vector3, Matrix4 } from 'three'
import Assembly from '../symmetry/assembly'
import Structure from '../structure/structure'
import Unitcell, { UnitcellParams } from '../symmetry/unitcell'

const reInteger = /^[1-9]$/
const reWhitespace = /\s+/
const bondTypes: {[k: string]: number} = {
  '1': 1,
  '2': 2,
  '3': 3,
  'am': 1, // amide
  'ar': 1, // aromatic
  'du': 1, // dummy
  'un': 1, // unknown
  'nc': 0 // not connected
}

interface IAccessStructuresMolData {
  mol2: string;
  spacegroupOperators: string;
  cellDimensions: string;
}

class CCDCMol2Parser extends StructureParser {
  get type () { return 'mol2' }

  _parse () {

    const jsonData = JSON.parse(this.streamer.asText()) as IAccessStructuresMolData


    // http://paulbourke.net/dataformats/mol2/

    if (Debug) Log.time('Mol2Parser._parse ' + this.name)

    const s = this.structure
    const sb = this.structureBuilder

    const firstModelOnly = this.firstModelOnly
    const asTrajectory = this.asTrajectory

    const frames = s.frames
    let doFrames = false
    let currentFrame: Float32Array, currentCoord: number

    const atomMap = s.atomMap
    const atomStore = s.atomStore
    atomStore.resize(Math.round(this.streamer.data.length / 60))
    atomStore.addField('partialCharge', 1, 'float32')

    let idx = 0
    let moleculeLineNo = 0
    let modelAtomIdxStart = 0
    let modelIdx = -1
    let numAtoms = 0

    let currentRecordType = 0
    let moleculeRecordType = 1
    let atomRecordType = 2
    let bondRecordType = 3
    let crysinRecordType = 4

    const unitcellDict: Partial<{
      origx: Matrix4
      scale: Matrix4
      a: number
      b: number
      c: number
      alpha: number
      beta: number
      gamma: number
      spacegroup: string
    }> = {}

    const ap1 = s.getAtomProxy()
    const ap2 = s.getAtomProxy()

    function _parseChunkOfLines (_i: number, _n: number, lines: string[]) {
      for (let i = _i; i < _n; ++i) {
        const line = lines[ i ].trim()

        if (line === '' || line[ 0 ] === '#') continue

        if (line[ 0 ] === '@') {
          if (line === '@<TRIPOS>MOLECULE') {
            currentRecordType = moleculeRecordType
            moleculeLineNo = 0

            ++modelIdx
          } else if (line === '@<TRIPOS>ATOM') {
            currentRecordType = atomRecordType
            modelAtomIdxStart = atomStore.count

            if (asTrajectory) {
              currentCoord = 0
              currentFrame = new Float32Array(numAtoms * 3)
              frames.push(currentFrame)

              if (modelIdx > 0) doFrames = true
            }
          } else if (line === '@<TRIPOS>BOND') {
            currentRecordType = bondRecordType
          } else if (line === '@<TRIPOS>CRYSIN') {
              currentRecordType = crysinRecordType
          } else {
            currentRecordType = 0
          }
        } else if (currentRecordType === moleculeRecordType) {
          if (moleculeLineNo === 0) {
            s.title = line
            s.id = line
          } else if (moleculeLineNo === 1) {
            const ls = line.split(reWhitespace)
            numAtoms = parseInt(ls[ 0 ])
            // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
          } else if (moleculeLineNo === 2) {

            // const molType = line;
            // SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE

          } else if (moleculeLineNo === 3) {

            // const chargeType = line;
            // NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL,
            // PULLMAN, GAUSS80_CHARGES, AMPAC_CHARGES,
            // MULLIKEN_CHARGES, DICT_ CHARGES, MMFF94_CHARGES,
            // USER_CHARGES

          } else if (moleculeLineNo === 4) {

            // const statusBits = line;

          } else if (moleculeLineNo === 5) {

            // const molComment = line;

          }

          ++moleculeLineNo
        } else if (currentRecordType === atomRecordType) {
          const ls = line.split(reWhitespace)

          if (firstModelOnly && modelIdx > 0) continue

          const x = parseFloat(ls[ 2 ])
          const y = parseFloat(ls[ 3 ])
          const z = parseFloat(ls[ 4 ])

          if (asTrajectory) {
            const j = currentCoord * 3

            currentFrame[ j + 0 ] = x
            currentFrame[ j + 1 ] = y
            currentFrame[ j + 2 ] = z

            currentCoord += 1

            if (doFrames) continue
          }

          const serial = ls[ 0 ]
          const atomname = ls[ 1 ]
          const element = ls[ 5 ].split('.')[ 0 ]
          const resno = ls[ 6 ] ? parseInt(ls[ 6 ]) : 1
          const resname = ls[ 7 ] ? ls[ 7 ] : ''
          const partialCharge = ls[ 8 ] ? parseFloat(ls[ 8 ]) : 0.0

          atomStore.growIfFull()
          atomStore.atomTypeId[ idx ] = atomMap.add(atomname, element)

          atomStore.x[ idx ] = x
          atomStore.y[ idx ] = y
          atomStore.z[ idx ] = z
          atomStore.serial[ idx ] = serial
          atomStore.partialCharge[ idx ] = partialCharge

          sb.addAtom(modelIdx, '', '', resname, resno, 1)

          idx += 1
        } else if (currentRecordType === bondRecordType) {
          if (firstModelOnly && modelIdx > 0) continue
          if (asTrajectory && modelIdx > 0) continue

          const ls = line.split(reWhitespace)

          // ls[ 0 ] is bond id
          ap1.index = parseInt(ls[ 1 ]) - 1 + modelAtomIdxStart
          ap2.index = parseInt(ls[ 2 ]) - 1 + modelAtomIdxStart
          const order = bondTypes[ ls[ 3 ] ]

          s.bondStore.addBond(ap1, ap2, order)
        }  else if (currentRecordType === crysinRecordType) {
          const ls = line.split(reWhitespace);
  
          const aLength = parseFloat(ls[0]);
          const bLength = parseFloat(ls[1]);
          const cLength = parseFloat(ls[2]);
         
          const alpha = parseFloat(ls[3]);
          const beta = parseFloat(ls[4]);
          const gamma = parseFloat(ls[5]);
  
          unitcellDict.a = aLength;
          unitcellDict.b = bLength;
          unitcellDict.c = cLength;
          unitcellDict.alpha = alpha;
          unitcellDict.beta = beta;
          unitcellDict.gamma = gamma;
        }
      }
    }

    const lines = jsonData.mol2.split('\n');
    _parseChunkOfLines(0, lines.length, lines)

    if (unitcellDict.a !== undefined) { 
      s.unitcell = new Unitcell(unitcellDict as UnitcellParams) 
    } else { 
      s.unitcell = undefined 
    } 

    const operators:string[][] = [];
    const spacegroupOperators = (spacegroupOperators: string) =>{
            // x,y,z;-x,-y,1/2+z;1/2-x,1/2+y,1/2+z;1/2+x,1/2-y,z
      spacegroupOperators.split(";").forEach(op => {
        const semop:string[] = [];
        op.split(",").forEach(o => {
          semop.push(o.toUpperCase());
        });
        operators.push(semop);
      });

      this._buildUnitcellAssembly(s, operators);
    }

    sb.finalize()
    s.finalizeAtoms()
    calculateChainnames(s)
    calculateBondsWithin(s, true)
    calculateBondsBetween(s, true)
    s.finalizeBonds()
    assignResidueTypeBonds(s)
    calculateSecondaryStructure(s)
    spacegroupOperators(jsonData.spacegroupOperators);

    if (Debug) Log.timeEnd('Mol2Parser._parse ' + this.name)
  }

  _setSymmetryOperations (operators: string[][]) {

    const matrixDict: { [k: string]: Matrix4 } = {}
  
    if (operators === undefined) {
      console.warn(`spacegroup '${operators}' not found`)
      return matrixDict
    }
  
    operators.forEach(function (symop) {
      let row = 0
      const matrix = new Matrix4().set(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1
      )
      const me = matrix.elements
  
      matrixDict[ symop.toString() ] = matrix
  
      symop.forEach(function (elm) {
        let negate = false
        let denominator = false
  
        for (let i = 0, n = elm.length; i < n; ++i) {
          const c = elm[ i ]
  
          if (c === '-') {
            negate = true
          } else if (c === '+') {
            negate = false
          } else if (c === '/') {
            denominator = true
          } else if (c === 'X') {
            me[ 0 + row ] = negate ? -1 : 1
          } else if (c === 'Y') {
            me[ 4 + row ] = negate ? -1 : 1
          } else if (c === 'Z') {
            me[ 8 + row ] = negate ? -1 : 1
          } else if (reInteger.test(c)) {
            const integer = parseInt(c)
            if (denominator) {
              me[ 12 + row ] /= integer
            } else {
              me[ 12 + row ] = integer
            }
          } else {
            Log.warn(`getSymmetryOperations: unknown token '${c}'`)
          }
        }
  
        row += 1
      })
    })
  
    return matrixDict
  }

  _buildUnitcellAssembly (structure: Structure, operators: string[][]) {
    if (!structure.unitcell) return
  
    if (Debug) Log.time('buildUnitcellAssembly')
  
    const uc = structure.unitcell
  
    const structureCenterFrac = structure.center.clone().applyMatrix4(uc.cartToFrac)
    const centerFrac = structureCenterFrac.clone().floor()
    const symopDict: { [K: string]: Matrix4 } = this._setSymmetryOperations(operators)
  
    const centerFracSymop = new Vector3()
    const positionFracSymop = new Vector3()
  
    function getMatrixList (shift?: Vector3) {
      const matrixList: Matrix4[] = []
  
      Object.keys(symopDict).forEach(function (name) {
        const m = symopDict[ name ].clone()
  
        centerFracSymop.copy(structureCenterFrac).applyMatrix4(m).floor()
        positionFracSymop.setFromMatrixPosition(m)
        positionFracSymop.sub(centerFracSymop)
        positionFracSymop.add(centerFrac)
  
        if (shift) positionFracSymop.add(shift)
  
        m.setPosition(positionFracSymop)
        m.multiplyMatrices(uc.fracToCart, m)
        m.multiply(uc.cartToFrac)
  
        matrixList.push(m)
      })
  
      return matrixList
    }
  
    const unitcellAssembly = new Assembly('UNITCELL')
    const unitcellMatrixList = getMatrixList()
    const ncsMatrixList: Matrix4[] = []
    if (structure.biomolDict.NCS) {
      ncsMatrixList.push(
        new Matrix4(), ...structure.biomolDict.NCS.partList[ 0 ].matrixList
      )
      const ncsUnitcellMatrixList: Matrix4[] = []
      unitcellMatrixList.forEach(sm => {
        ncsMatrixList.forEach(nm => {
          ncsUnitcellMatrixList.push(sm.clone().multiply(nm))
        })
      })
      unitcellAssembly.addPart(ncsUnitcellMatrixList)
    } else {
      unitcellAssembly.addPart(unitcellMatrixList)
    }
  
    const vec = new Vector3()
    const supercellAssembly = new Assembly('SUPERCELL')
    const supercellMatrixList = Array.prototype.concat.call(
      getMatrixList(vec.set(1, 0, 0)),  // 655
      getMatrixList(vec.set(0, 1, 0)),  // 565
      getMatrixList(vec.set(0, 0, 1)),  // 556
  
      getMatrixList(vec.set(-1, 0, 0)),  // 455
      getMatrixList(vec.set(0, -1, 0)),  // 545
      getMatrixList(vec.set(0, 0, -1)),  // 554
  
      getMatrixList(vec.set(1, 1, 0)),  // 665
      getMatrixList(vec.set(1, 0, 1)),  // 656
      getMatrixList(vec.set(0, 1, 1)),  // 566
  
      getMatrixList(vec.set(-1, -1, 0)),  // 445
      getMatrixList(vec.set(-1, 0, -1)),  // 454
      getMatrixList(vec.set(0, -1, -1)),  // 544
  
      getMatrixList(vec.set(1, -1, -1)),  // 644
      getMatrixList(vec.set(1, 1, -1)),  // 664
      getMatrixList(vec.set(1, -1, 1)),  // 646
      getMatrixList(vec.set(-1, 1, 1)),  // 466
      getMatrixList(vec.set(-1, -1, 1)),  // 446
      getMatrixList(vec.set(-1, 1, -1)),  // 464
  
      getMatrixList(vec.set(0, 1, -1)),  // 564
      getMatrixList(vec.set(0, -1, 1)),  // 546
      getMatrixList(vec.set(1, 0, -1)),  // 654
      getMatrixList(vec.set(-1, 0, 1)),  // 456
      getMatrixList(vec.set(1, -1, 0)),  // 645
      getMatrixList(vec.set(-1, 1, 0)),  // 465
  
      getMatrixList(),  // 555
      getMatrixList(vec.set(1, 1, 1)),  // 666
      getMatrixList(vec.set(-1, -1, -1))   // 444
    )
    if (structure.biomolDict.NCS) {
      const ncsSupercellMatrixList: Matrix4[] = []
      supercellMatrixList.forEach(function (sm: Matrix4) {
        ncsMatrixList.forEach(function (nm) {
          ncsSupercellMatrixList.push(sm.clone().multiply(nm))
        })
      })
      supercellAssembly.addPart(ncsSupercellMatrixList)
    } else {
      supercellAssembly.addPart(supercellMatrixList)
    }
  
    structure.biomolDict.UNITCELL = unitcellAssembly
    structure.biomolDict.SUPERCELL = supercellAssembly
  
    if (Debug) Log.timeEnd('buildUnitcellAssembly')
  }
}

ParserRegistry.add('ccdcmol2', CCDCMol2Parser)

export default CCDCMol2Parser



