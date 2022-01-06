/**
 * @file Mol2 Parser
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @private
 */

import { Matrix4 } from 'three'
import { Debug, Log, ParserRegistry } from '../globals'
import {
  assignResidueTypeBonds,
  calculateChainnames, calculateSecondaryStructure,
  calculateBondsBetween, calculateBondsWithin, buildUnitcellAssembly
} from '../structure/structure-utils'
import Unitcell, { UnitcellParams } from '../symmetry/unitcell'
import StructureParser from './structure-parser'

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

const spaceGroupLookup: string[][] = [
  ["P 1", "", "", "", "", ""],
  ["P -1", "", "", "", "", ""],
  ["P 2", "", "P 1 1 2", "", "", ""],
  ["P 21", "", "P 1 1 21", "", "", ""],
  ["C 2", "A 2", "B 1 1 2", "A 1 1 2", "", ""],
  ["P m", "", "P 1 1 m", "", "", ""],
  ["P c", "P a", "P 1 1 b", "P 1 1 a", "", ""],
  ["C m", "A m", "B 1 1 m", "A 1 1 m", "", ""],
  ["C c", "A a", "B 1 1 b", "A 1 1 a", "", ""],
  ["P 2/m", "", "P 1 1 2/m", "", "", ""],
  ["P 21/m", "", "P 1 1 21/m", "", "", ""],
  ["C 2/m", "A 2/m", "B 1 1 2/m", "A 1 1 2/m", "", ""],
  ["P 2/c", "P 2/a", "P 1 1 2/b", "P 1 1 2/a", "", ""],
  ["P 21/c", "P 21/a", "P 1 1 21/b", "P 1 1 21/a", "", ""],
  ["C 2/c", "A 2/a", "B 1 1 2/b", "A 1 1 2/a", "", ""],
  ["P 2 2 2", "", "", "", "", ""],
  ["P 2 2 21", "P 21 2 2", "P 2 21 2", "", "", ""],
  ["P 21 21 2", "P 2 21 21", "P 21 2 21", "", "", ""],
  ["P 21 21 21", "", "", "", "", ""],
  ["C 2 2 21", "A 21 2 2", "B 2 21 2", "", "", ""],
  ["C 2 2 2", "A 2 2 2", "B 2 2 2", "", "", ""],
  ["F 2 2 2", "", "", "", "", ""],
  ["I 2 2 2", "", "", "", "", ""],
  ["I 21 21 21", "", "", "", "", ""],
  ["P m m 2", "P 2 m m", "P m 2 m", "", "", ""],
  ["P m c 21", "P 21 m a", "P b 21 m", "P m 21 b", "P c m 21", "P 21 a m"],
  ["P c c 2", "P 2 a a", "P b 2 b", "", "", ""],
  ["P m a 2", "P 2 m b", "P c 2 m", "P m 2 a", "P b m 2", "P 2 c m"],
  ["P c a 21", "P 21 a b", "P c 21 b", "P b 21 a", "P b c 21", "P 21 c a"],
  ["P n c 2", "P 2 n a", "P b 2 n", "P n 2 b", "P c n 2", "P 2 a n"],
  ["P m n 21", "P 21 m n", "P n 21 m", "P m 21 n", "P n m 21", "P 21 n m"],
  ["P b a 2", "P 2 c b", "P c 2 a", "", "", ""],
  ["P n a 21", "P 21 n b", "P c 21 n", "P n 21 a", "P b n 21", "P 21 c n"],
  ["P n n 2", "P 2 n n", "P n 2 n", "", "", ""],
  ["C m m 2", "A 2 m m", "B m 2 m", "", "", ""],
  ["C m c 21", "A 21 m a", "B b 21 m", "B m 21 b", "C c m 21", "A 21 a m"],
  ["C c c 2", "A 2 a a", "B b 2 b", "", "", ""],
  ["A m m 2", "B 2 m m", "C m 2 m", "A m 2 m", "B m m 2", "C 2 m m"],
  ["A b m 2", "B 2 c m", "C m 2 a", "A c 2 m", "B m a 2", "C 2 m b"],
  ["A m a 2", "B 2 m b", "C c 2 m", "A m 2 a", "B b m 2", "C 2 c m"],
  ["A b a 2", "B 2 c b", "C c 2 a", "A c 2 a", "B b a 2", "C 2 c b"],
  ["F m m 2", "F 2 m m", "F m 2 m", "", "", ""],
  ["F d d 2", "F 2 d d", "F d 2 d", "", "", ""],
  ["I m m 2", "I 2 m m", "I m 2 m", "", "", ""],
  ["I b a 2", "I 2 c b", "I c 2 a", "", "", ""],
  ["I m a 2", "I 2 m b", "I c 2 m", "I m 2 a", "I b m 2", "I 2 c m"],
  ["P m m m", "", "", "", "", ""],
  ["P n n n", "", "", "", "", ""],
  ["P c c m", "P m a a", "P b m b", "", "", ""],
  ["P b a n", "P n c b", "P c n a", "", "", ""],
  ["P m m a", "P b m m", "P m c m", "P m a m", "P m m b", "P c m m"],
  ["P n n a", "P b n n", "P n c n", "P n a n", "P n n b", "P c n n"],
  ["P m n a", "P b m n", "P n c m", "P m a n", "P n m b", "P c n m"],
  ["P c c a", "P b a a", "P b c b", "P b a b", "P c c b", "P c a a"],
  ["P b a m", "P m c b", "P c m a", "", "", ""],
  ["P c c n", "P n a a", "P b n b", "", "", ""],
  ["P b c m", "P m c a", "P b m a", "P c m b", "P c a m", "P m a b"],
  ["P n n m", "P m n n", "P n m n", "", "", ""],
  ["P m m n", "P n m m", "P m n m", "", "", ""],
  ["P b c n", "P n c a", "P b n a", "P c n b", "P c a n", "P n a b"],
  ["P b c a", "", "P c a b", "", "", ""],
  ["P n m a", "P b n m", "P m c n", "P n a m", "P m n b", "P c m n"],
  ["C m c m", "A m m a", "B b m m", "B m m b", "C c m m", "A m a m"],
  ["C m c a", "A b m a", "B b c m", "B m a b", "C c m b", "A c a m"],
  ["C m m m", "A m m m", "B m m m", "", "", ""],
  ["C c c m", "A m a a", "B b m b", "", "", ""],
  ["C m m a", "A b m m", "B m c m", "B m a m", "C m m b", "A c m m"],
  ["C c c a", "A b a a", "B b c b", "B b a b", "C c c b", "A c a a"],
  ["F m m m", "", "", "", "", ""],
  ["F d d d", "", "", "", "", ""],
  ["I m m m", "", "", "", "", ""],
  ["I b a m", "I m c b", "I c m a", "", "", ""],
  ["I b c a", "", "I c a b", "", "", ""],
  ["I m m a", "I b m m", "I m c m", "I m a m", "I m m b", "I c m m"],
  ["P 4", "C 4", "", "", "", ""],
  ["P 41", "C 41", "", "", "", ""],
  ["P 42", "C 42", "", "", "", ""],
  ["P 43", "C 43", "", "", "", ""],
  ["I 4", "F 4", "", "", "", ""],
  ["I 41", "F 41", "", "", "", ""],
  ["P -4", "C -4", "", "", "", ""],
  ["I -4", "F -4", "", "", "", ""],
  ["P 4/m", "C 4/m", "", "", "", ""],
  ["P 42/m", "C 42/m", "", "", "", ""],
  ["P 4/n", "C 4/a", "", "", "", ""],
  ["P 42/n", "C 42/a", "", "", "", ""],
  ["I 4/m", "F 4/m", "", "", "", ""],
  ["I 41/a", "F 41/d", "", "", "", ""],
  ["P 4 2 2", "C 4 2 2", "", "", "", ""],
  ["P 4 21 2", "C 4 2 21", "", "", "", ""],
  ["P 41 2 2", "C 41 2 2", "", "", "", ""],
  ["P 41 21 2", "C 41 2 21", "", "", "", ""],
  ["P 42 2 2", "C 42 2 2", "", "", "", ""],
  ["P 42 21 2", "C 42 2 21", "", "", "", ""],
  ["P 43 2 2", "C 43 2 2", "", "", "", ""],
  ["P 43 21 2", "C 43 2 21", "", "", "", ""],
  ["I 4 2 2", "F 4 2 2", "", "", "", ""],
  ["I 41 2 2", "F 41 2 2", "", "", "", ""],
  ["P 4 m m", "C 4 m m", "", "", "", ""],
  ["P 4 b m", "C 4 m b", "", "", "", ""],
  ["P 42 c m", "C 42 m c", "", "", "", ""],
  ["P 42 n m", "C 42 m n", "", "", "", ""],
  ["P 4 c c", "C 4 c c", "", "", "", ""],
  ["P 4 n c", "C 4 c n", "", "", "", ""],
  ["P 42 m c", "C 42 c m", "", "", "", ""],
  ["P 42 b c", "C 42 c b", "", "", "", ""],
  ["I 4 m m", "F 4 m m", "", "", "", ""],
  ["I 4 c m", "F 4 m c", "", "", "", ""],
  ["I 41 m d", "F 41 d m", "", "", "", ""],
  ["I 41 c d", "F 41 d c", "", "", "", ""],
  ["P -4 2 m", "C -4 m 2", "", "", "", ""],
  ["P -4 2 c", "C -4 c 2", "", "", "", ""],
  ["P -4 21 m", "C -4 m 21", "", "", "", ""],
  ["P -4 21 c", "C -4 c 21", "", "", "", ""],
  ["P -4 m 2", "C -4 2 m", "", "", "", ""],
  ["P -4 c 2", "C -4 2 c", "", "", "", ""],
  ["P -4 b 2", "C -4 2 b", "", "", "", ""],
  ["P -4 n 2", "C -4 2 n", "", "", "", ""],
  ["I -4 m 2", "F -4 2 m", "", "", "", ""],
  ["I -4 c 2", "F -4 2 c", "", "", "", ""],
  ["I -4 2 m", "F -4 m 2", "", "", "", ""],
  ["I -4 2 d", "F -4 d 2", "", "", "", ""],
  ["P 4/m m m", "C 4/m m m", "", "", "", ""],
  ["P 4/m c c", "C 4/m c c", "", "", "", ""],
  ["P 4/n b m", "C 4/a m b", "", "", "", ""],
  ["P 4/n n c", "C 4/a c n", "", "", "", ""],
  ["P 4/m b m", "C 4/m m b", "", "", "", ""],
  ["P 4/m n c", "C 4/m c n", "", "", "", ""],
  ["P 4/n m m", "C 4/a m m", "", "", "", ""],
  ["P 4/n c c", "C 4/a c c", "", "", "", ""],
  ["P 42/m m c", "C 42/m c m", "", "", "", ""],
  ["P 42/m c m", "C 42/m m c", "", "", "", ""],
  ["P 42/n b c", "C 42/a c b", "", "", "", ""],
  ["P 42/n n m", "C 42/a m n", "", "", "", ""],
  ["P 42/m b c", "C 42/m c b", "", "", "", ""],
  ["P 42/m n m", "C 42/m m n", "", "", "", ""],
  ["P 42/n m c", "C 42/a c m", "", "", "", ""],
  ["P 42/n c m", "C 42/a m c", "", "", "", ""],
  ["I 4/m m m", "F 4/m m m", "", "", "", ""],
  ["I 4/m c m", "F 4/m m c", "", "", "", ""],
  ["I 41/a m d", "F 41/d d m", "", "", "", ""],
  ["I 41/a c d", "F 41/d d c", "", "", "", ""],
  ["P 3", "", "", "", "", ""],
  ["P 31", "", "", "", "", ""],
  ["P 32", "", "", "", "", ""],
  ["R 3", "", "", "", "", ""],
  ["P -3", "", "", "", "", ""],
  ["R -3", "", "", "", "", ""],
  ["P 3 1 2", "", "", "", "", ""],
  ["P 3 2 1", "", "", "", "", ""],
  ["P 31 1 2", "", "", "", "", ""],
  ["P 31 2 1", "", "", "", "", ""],
  ["P 32 1 2", "", "", "", "", ""],
  ["P 32 2 1", "", "", "", "", ""],
  ["R 3 2", "", "", "", "", ""],
  ["P 3 m 1", "", "", "", "", ""],
  ["P 3 1 m", "", "", "", "", ""],
  ["P 3 c 1", "", "", "", "", ""],
  ["P 3 1 c", "", "", "", "", ""],
  ["R 3 m", "", "", "", "", ""],
  ["R 3 c", "", "", "", "", ""],
  ["P -3 1 m", "", "", "", "", ""],
  ["P -3 1 c", "", "", "", "", ""],
  ["P -3 m 1", "", "", "", "", ""],
  ["P -3 c 1", "", "", "", "", ""],
  ["R -3 m", "", "", "", "", ""],
  ["R -3 c", "", "", "", "", ""],
  ["P 6", "", "", "", "", ""],
  ["P 61", "", "", "", "", ""],
  ["P 65", "", "", "", "", ""],
  ["P 62", "", "", "", "", ""],
  ["P 64", "", "", "", "", ""],
  ["P 63", "", "", "", "", ""],
  ["P -6", "", "", "", "", ""],
  ["P 6/m", "", "", "", "", ""],
  ["P 63/m", "", "", "", "", ""],
  ["P 6 2 2", "", "", "", "", ""],
  ["P 61 2 2", "", "", "", "", ""],
  ["P 65 2 2", "", "", "", "", ""],
  ["P 62 2 2", "", "", "", "", ""],
  ["P 64 2 2", "", "", "", "", ""],
  ["P 63 2 2", "", "", "", "", ""],
  ["P 6 m m", "", "", "", "", ""],
  ["P 6 c c", "", "", "", "", ""],
  ["P 63 c m", "", "", "", "", ""],
  ["P 63 m c", "", "", "", "", ""],
  ["P -6 m 2", "", "", "", "", ""],
  ["P -6 c 2", "", "", "", "", ""],
  ["P -6 2 m", "", "", "", "", ""],
  ["P -6 2 c", "", "", "", "", ""],
  ["P 6/m m m", "", "", "", "", ""],
  ["P 6/m c c", "", "", "", "", ""],
  ["P 63/m c m", "", "", "", "", ""],
  ["P 63/m m c", "", "", "", "", ""],
  ["P 2 3", "", "", "", "", ""],
  ["F 2 3", "", "", "", "", ""],
  ["I 2 3", "", "", "", "", ""],
  ["P 21 3", "", "", "", "", ""],
  ["I 21 3", "", "", "", "", ""],
  ["P m -3", "", "", "", "", ""],
  ["P n -3", "", "", "", "", ""],
  ["F m -3", "", "", "", "", ""],
  ["F d -3", "", "", "", "", ""],
  ["I m -3", "", "", "", "", ""],
  ["P a -3", "", "", "", "", ""],
  ["I a -3", "", "", "", "", ""],
  ["P 4 3 2", "", "", "", "", ""],
  ["P 42 3 2", "", "", "", "", ""],
  ["F 4 3 2", "", "", "", "", ""],
  ["F 41 3 2", "", "", "", "", ""],
  ["I 4 3 2", "", "", "", "", ""],
  ["P 43 3 2", "", "", "", "", ""],
  ["P 41 3 2", "", "", "", "", ""],
  ["I 41 3 2", "", "", "", "", ""],
  ["P -4 3 m", "", "", "", "", ""],
  ["F -4 3 m", "", "", "", "", ""],
  ["I -4 3 m", "", "", "", "", ""],
  ["P -4 3 n", "", "", "", "", ""],
  ["F -4 3 c", "", "", "", "", ""],
  ["I -4 3 d", "", "", "", "", ""],
  ["P m -3 m", "", "", "", "", ""],
  ["P n -3 n", "", "", "", "", ""],
  ["P m -3 n", "", "", "", "", ""],
  ["P n -3 m", "", "", "", "", ""],
  ["F m -3 m", "", "", "", "", ""],
  ["F m -3 c", "", "", "", "", ""],
  ["F d -3 m", "", "", "", "", ""],
  ["F d -3 c", "", "", "", "", ""],
  ["I m -3 m", "", "", "", "", ""],
  ["I a -3 d", "", "", "", "", ""],
];

class Mol2Parser extends StructureParser {
  get type () { return 'mol2' }

  _parse () {
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

    const ap1 = s.getAtomProxy()
    const ap2 = s.getAtomProxy()

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
          const spaceGroup = parseFloat(ls[6]);
          const setting = parseFloat(ls[7]);

          unitcellDict.a = aLength;
          unitcellDict.b = bLength;
          unitcellDict.c = cLength;
          unitcellDict.alpha = alpha;
          unitcellDict.beta = beta;
          unitcellDict.gamma = gamma;
        
          unitcellDict.spacegroup = lookupSpaceGroup(spaceGroup, setting);
        }
      }
    }

    this.streamer.eachChunkOfLines(function (lines/*, chunkNo, chunkCount */) {
      _parseChunkOfLines(0, lines.length, lines)
    })

    if (unitcellDict.a !== undefined) { 
      s.unitcell = new Unitcell(unitcellDict as UnitcellParams) 
    } else { 
      s.unitcell = undefined 
    } 
    
    sb.finalize()
    s.finalizeAtoms()
    calculateChainnames(s)
    calculateBondsWithin(s, true)
    calculateBondsBetween(s, true)
    s.finalizeBonds()
    assignResidueTypeBonds(s)
    calculateSecondaryStructure(s)
    buildUnitcellAssembly(s)

    if (Debug) Log.timeEnd('Mol2Parser._parse ' + this.name)
  }
}

ParserRegistry.add('mol2', Mol2Parser)

export default Mol2Parser
function lookupSpaceGroup(spaceGroup: number, setting: number): string | undefined {
    if(spaceGroup > 230 && setting == 1){
      switch(spaceGroup){
        case 231:
          spaceGroup = 14;
          setting = 2;
          break;
        case 232:
          spaceGroup = 146;
          break;
        case 233:
          spaceGroup = 148;
        break;
        case 234:
          spaceGroup = 155;
          break;
        case 235:
          spaceGroup = 160;
          break;
        case 236:
          spaceGroup = 161;
          break;
        case 237:
          spaceGroup = 166;
          break;
        case 238:
          spaceGroup = 167;
          break;
        default:
          break;
      }
  }
  return spaceGroupLookup[spaceGroup - 1][setting - 1];
}