import Mol2Parser from "./mol2-parser";
import { Debug, Log, ParserRegistry } from '../globals'
import { Vector3, Matrix4 } from 'three'
import Assembly from '../symmetry/assembly'
import Structure from '../structure/structure'
import Unitcell, { UnitcellParams } from '../symmetry/unitcell'
import StringStreamer from "../streamer/string-streamer";
import Streamer from "../streamer/streamer";
import { StructureParserParameters } from "./structure-parser";

const reInteger = /^[1-9]$/
const reWhitespace = /\s+/

interface IAccessStructuresMolData {
    mol2: string;
    spacegroupOperators: string;
    cellDimensions: string;
}

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
  
class CCDCMol2Parser {
  streamer: Streamer
  name: string
  path: string
  [k: string]: any

    get type () { return 'ccdcmol2' }
    constructor(streamer: Streamer, params?: Partial<StructureParserParameters>) {      
        this.streamer = streamer;
        this.params = params;     
    }

    parse () {
      return this.streamer.read().then(() => {
        const data = JSON.parse(this.streamer.asText()) as IAccessStructuresMolData;
        const mol2Parser = new Mol2Parser(new StringStreamer(data.mol2), this.params);
        return mol2Parser.parse().then((structure)=> {
          buildUnitCell(structure, data.cellDimensions);
          spacegroupOperators(structure, data.spacegroupOperators);
          return structure
        });
      })
    }
}

ParserRegistry.add('ccdcmol2', CCDCMol2Parser)

export default CCDCMol2Parser


/**
 * Builds unit cell based on the cell dimensions in the file
 * @param s 
 * @param cellDimensions 
 */
function buildUnitCell(s: Structure, cellDimensions: string): void {
  if(cellDimensions)
  {
    const ls =  cellDimensions.split(reWhitespace);

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

  if (unitcellDict.a !== undefined) { 
    s.unitcell = new Unitcell(unitcellDict as UnitcellParams) 
  } else { 
    s.unitcell = undefined 
  }
}

/**
 * Split operators into a 2D array in preparation for generating the unit cell assembly
 * @param spacegroupOperators 
 */    
function spacegroupOperators(s: Structure, spacegroupOperators: string): void {
  const operators: string[][] = [];
  spacegroupOperators.split(";").forEach(op => {
    operators.push(op.split(",").map(v => v.toUpperCase()));
  });
  
  buildUnitcellAssembly(s, operators);
}


function setSymmetryOperations (operators: string[][]): { [k: string]: Matrix4 } {

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

function buildUnitcellAssembly(structure: Structure, operators: string[][]): void {
  if (!structure.unitcell) return
  if (Debug) Log.time('buildUnitcellAssembly')
  
  const uc = structure.unitcell

  const structureCenterFrac = structure.center.clone().applyMatrix4(uc.cartToFrac)
  const centerFrac = structureCenterFrac.clone().floor()
  const symopDict: { [K: string]: Matrix4 } = setSymmetryOperations(operators)

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