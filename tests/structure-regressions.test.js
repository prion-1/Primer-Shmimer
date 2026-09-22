'use strict';

const {
    test,
    before,
    afterEach
} = require('node:test');
const assert = require('node:assert/strict');

global.window = globalThis;
global.ThalModule = require('../primer3_wasm/dist/thal.js');

require('../thal_bridge.js');
require('../thermodynamics.js');
require('../algorithms.js');

const BioAlgorithms = global.BioAlgorithms;
const structureOptions = Object.freeze({
    naConc: 0.05,
    mgConc: 0.0015,
    dNTPConc: 0.0008,
    kConc: 0,
    trisConc: 0,
    primerConc: 200e-9,
    evalTempK: 333.15,
    maxLoopSize: 30
});

const neutralThermo = Object.freeze({ tm: 60 });
const neutralMetrics = Object.freeze({
    end3StabilityDG: 0,
    gcClampCount: 2,
    hasRuns: false
});
const neutralPrimer = 'ATGCGATCGATCGATCGCAT';
const emptyStructureResult = Object.freeze({
    dG: 0,
    dH: 0,
    dS: 0,
    Tm: 0,
    alignEnd1: -1,
    alignEnd2: -1,
    structure: null,
    structureFound: false,
    message: ''
});

let realPrimer3;

before(async () => {
    const ready = await global.Primer3.init();
    assert.equal(ready, true, 'Primer3 WASM should initialize for regression tests');
    realPrimer3 = global.Primer3;
});

afterEach(() => {
    global.Primer3 = realPrimer3;
});

test('structure analyzers are part of the public BioAlgorithms API', () => {
    assert.equal(typeof BioAlgorithms.analyzeHairpinStructure, 'function');
    assert.equal(typeof BioAlgorithms.analyzeDimerStructure, 'function');
});

test('the screenshot GNA mini-hairpin keeps its measured high Tm and structural annotation', () => {
    const result = BioAlgorithms.findBestHairpin('GAAGAGCGGAGCC', structureOptions);

    assert.equal(result.structureFound, true);
    assert.ok(Math.abs(result.Tm - 69.01543273233278) < 1e-9);
    assert.ok(Math.abs(result.dG - (-0.38204845409117477)) < 1e-9);
    assert.equal(result.basePairCount, 2);
    assert.equal(result.specialTriloop, true);
    assert.equal(result.isThreePrimeExtensible, false);
});

test('hairpin pair count counts pair tokens rather than the span across a bulge', () => {
    const result = BioAlgorithms.findBestHairpin('TCCAGCCTCGGCGCAAGGCC', structureOptions);

    assert.equal(result.structureFound, true);
    assert.match(result.structure, /\/-\//);
    assert.equal(result.basePairCount, 2);
    assert.equal(result.isThreePrimeExtensible, false);
});

test('hairpin extensibility requires the terminal 3-prime base itself to be paired', () => {
    const pairedEnd = BioAlgorithms.findBestHairpin('AAGCGGAGC', structureOptions);
    const unpairedTail = BioAlgorithms.findBestHairpin('AGCGGAGCA', structureOptions);

    assert.equal(pairedEnd.structureFound, true);
    assert.equal(pairedEnd.threePrimeRunLength, 2);
    assert.equal(pairedEnd.fivePrimeTemplateOverhang, 2);
    assert.equal(pairedEnd.isThreePrimeExtensible, true);

    assert.equal(unpairedTail.structureFound, true);
    assert.equal(unpairedTail.threePrimeRunLength, 0);
    assert.equal(unpairedTail.isThreePrimeExtensible, false);
});

test('dimer END requires adjacent terminal canonical pairs and a template overhang', () => {
    const positive = BioAlgorithms.findBestDimer('AAAAAGCGC', 'TTTTTGCGC', structureOptions);
    assert.equal(positive.END.structureFound, true);
    assert.equal(positive.END.basePairCount, 4);
    assert.equal(positive.END.isThreePrimeExtensible, true);

    const negative = BioAlgorithms.findBestDimer('GCGCGAAAA', 'TTTTGCGCG', structureOptions);
    assert.equal(negative.END.structureFound, false);
    assert.equal(negative.END.isThreePrimeExtensible, false);
    assert.match(negative.END.message, /No 3' extensible dimer found/);
});

test('cross-dimer END checks both primer orientations', () => {
    const calls = [];
    const reverseOnlyEnd = {
        dG: -6,
        dH: -30000,
        dS: -90,
        Tm: 40,
        alignEnd1: 2,
        alignEnd2: 2,
        structure: 'SEQ\tCG--\nSTR\tGCAA\n',
        structureFound: true,
        message: ''
    };

    global.Primer3 = {
        isReady: () => true,
        calcHairpin: () => emptyStructureResult,
        calcDimer: (seq1, seq2, options, mode) => {
            calls.push({ seq1, seq2, mode });
            if (mode === 'END' && seq1 === 'CG' && seq2 === 'AAAA') {
                return { ...reverseOnlyEnd };
            }
            return { ...emptyStructureResult };
        }
    };

    const result = BioAlgorithms.findBestDimer('AAAA', 'CG', structureOptions);
    const endCalls = calls.filter(call => call.mode === 'END');

    const evaluatedOrientations = new Set(
        endCalls.map(call => `${call.seq1}>${call.seq2}`)
    );
    assert.equal(evaluatedOrientations.has('AAAA>CG'), true);
    assert.equal(evaluatedOrientations.has('CG>AAAA'), true);
    assert.equal(result.END.structureFound, true);
    assert.equal(result.END.basePairCount, 2);
    assert.equal(result.END.isThreePrimeExtensible, true);
});

test('a terminal mismatch does not make a dimer extensible', () => {
    const mismatchEnd = {
        dG: -6,
        dH: -30000,
        dS: -90,
        Tm: 40,
        alignEnd1: 2,
        alignEnd2: 2,
        structure: 'SEQ\tCG--\nSTR\tGAAA\n',
        structureFound: true,
        message: ''
    };

    global.Primer3 = {
        isReady: () => true,
        calcHairpin: () => emptyStructureResult,
        calcDimer: (seq1, seq2, options, mode) => (
            mode === 'END' ? { ...mismatchEnd } : { ...emptyStructureResult }
        )
    };

    const result = BioAlgorithms.findBestDimer('CG', 'AAAA', structureOptions);

    assert.equal(result.END.structureFound, false);
    assert.equal(result.END.isThreePrimeExtensible, false);
    assert.match(result.END.message, /No 3' extensible dimer found/);
});

test('an unpaired 3-prime hairpin tail does not receive the self-priming penalty', () => {
    const baseHairpin = {
        dG: -2,
        dH: -20000,
        dS: -60,
        Tm: 50,
        structureFound: true,
        basePairCount: 3,
        specialTriloop: false
    };
    const dimer = {
        ANY: { ...emptyStructureResult },
        END: { ...emptyStructureResult, isThreePrimeExtensible: false }
    };
    const pairedEnd = {
        ...baseHairpin,
        structure: 'SEQ\t---///---\\\\\\\nSTR\tAAAGCGAAACGC\n',
        isThreePrimeExtensible: true
    };
    const unpairedTail = {
        ...baseHairpin,
        structure: 'SEQ\t---///---\\\\\\-\nSTR\tAAAGCGAAACGCA\n',
        isThreePrimeExtensible: false
    };

    const pairedReport = BioAlgorithms.getPrimerScoreDetails(
        neutralThermo,
        pairedEnd,
        dimer,
        neutralMetrics,
        neutralPrimer,
        60
    );
    const tailReport = BioAlgorithms.getPrimerScoreDetails(
        neutralThermo,
        unpairedTail,
        dimer,
        neutralMetrics,
        neutralPrimer,
        60
    );

    assert.equal(
        pairedReport.penalties.some(penalty => penalty.label === "Hairpin 3' self-priming"),
        true
    );
    assert.equal(
        tailReport.penalties.some(penalty => penalty.label === "Hairpin 3' self-priming"),
        false
    );
    assert.equal(pairedReport.score, 75);
    assert.equal(tailReport.score, 90);
});

test('dimer score deduplication follows sameAsGlobal rather than equal energies', () => {
    const hairpin = { ...emptyStructureResult, specialTriloop: false };
    const globalDimer = {
        ...emptyStructureResult,
        dG: -6,
        structure: 'global alignment',
        structureFound: true
    };
    const differentEnd = {
        ...emptyStructureResult,
        dG: -6,
        structure: 'different end alignment',
        structureFound: true,
        isThreePrimeExtensible: true,
        sameAsGlobal: false
    };
    const sameEnd = {
        ...globalDimer,
        isThreePrimeExtensible: true,
        sameAsGlobal: true
    };

    const distinctReport = BioAlgorithms.getPrimerScoreDetails(
        neutralThermo,
        hairpin,
        { ANY: globalDimer, END: differentEnd },
        neutralMetrics,
        neutralPrimer,
        60
    );
    const dedupedReport = BioAlgorithms.getPrimerScoreDetails(
        neutralThermo,
        hairpin,
        { ANY: globalDimer, END: sameEnd },
        neutralMetrics,
        neutralPrimer,
        60
    );

    assert.deepEqual(
        distinctReport.penalties
            .filter(penalty => penalty.label.includes('dimer'))
            .map(penalty => penalty.points),
        [10, 15]
    );
    assert.deepEqual(
        dedupedReport.penalties
            .filter(penalty => penalty.label.includes('dimer'))
            .map(penalty => penalty.points),
        [15]
    );
    assert.equal(distinctReport.score, 75);
    assert.equal(dedupedReport.score, 85);
});
