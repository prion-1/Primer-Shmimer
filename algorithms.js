// algorithms.js
(function() {
    const DB = window.ThermoDB;

    function comp(base) {
        const pairs = { A: 'T', T: 'A', C: 'G', G: 'C' };
        return pairs[base] || '';
    }

    function isPalindrome(seq) {
        for (let i = 0; i < seq.length / 2; i++) {
            if (seq[i] !== comp(seq[seq.length - 1 - i])) return false;
        }
        return true;
    }

    function calculateSaltCorrection(seqLen, fGC, naConc, kConc, trisConc, mgConc, dNTPConc, useOwczarzy = true) {
        const monovalent = naConc + kConc + trisConc / 2;
        // Guard: use 1e-11 instead of 0 when dNTP >= Mg2+ to avoid log(0)
        // in the Owczarzy lnMg term. Effectively zero free Mg2+.
        const freeMg = mgConc > dNTPConc ? (mgConc - dNTPConc) : 1e-11;

        let saltCorrEntropy = 0;
        let tmCorrOwczarzy = 0;

        if (!useOwczarzy) {
            const naEq = monovalent + 120 * Math.sqrt(freeMg);
            if (naEq > 0) saltCorrEntropy = 0.368 * (seqLen - 1) * Math.log(naEq);
            return { dS_corr: saltCorrEntropy, tm_inv_corr: null, R: -1 };
        }

        const divMonovRatio = monovalent === 0 ? 6.0 : Math.sqrt(freeMg) / monovalent;

        if (divMonovRatio < 0.22) {
            const lnMono = Math.log(monovalent);
            tmCorrOwczarzy =
                (((4.29 * fGC) - 3.95) * 1e-5 * lnMono) +
                (9.40e-6 * lnMono * lnMono);
        } else {
            let a = 3.92e-5;
            const b = -9.11e-6;
            const c = 6.26e-5;
            let d = 1.42e-5;
            const e = -4.82e-4;
            const f = 5.25e-4;
            let g = 8.31e-5;
            const lnMg = Math.log(freeMg);

            if (divMonovRatio < 6.0 && monovalent > 0) {
                const lnMono = Math.log(monovalent);
                a = 3.92e-5 * (0.843 - (0.352 * Math.sqrt(monovalent) * lnMono));
                d = 1.42e-5 * (1.279 - 4.03e-3 * lnMono - 8.03e-3 * lnMono * lnMono);
                g = 8.31e-5 * (0.486 - 0.258 * lnMono + 5.25e-3 * lnMono * lnMono * lnMono);
            }

            tmCorrOwczarzy =
                a +
                (b * lnMg) +
                fGC * (c + (d * lnMg));

            // Primer3's oligotm.c omits the (1/(2*(len-1)))*(e + f*lnMg + g*lnMg²)
            // term for primer-length sequences. Verified empirically: including it
            // degrades Tm agreement from MAE 0.026°C to 2.064°C against primer3-py.
        }

        return { dS_corr: 0, tm_inv_corr: tmCorrOwczarzy, R: divMonovRatio };
    }

    function calculateThermodynamics(seq, options = {}) {
        seq = seq.toUpperCase();
        if (seq.length < 2) return { tm: 0, dG: 0, dH: 0, dS: 0 };

        const naConc = options.naConc !== undefined ? options.naConc : 0.05;
        const mgConc = options.mgConc !== undefined ? options.mgConc : 0.0015;
        const dNTPConc = options.dNTPConc !== undefined ? options.dNTPConc : 0.0008;
        const kConc = options.kConc !== undefined ? options.kConc : 0;
        const trisConc = options.trisConc !== undefined ? options.trisConc : 0;
        const primerConc = options.primerConc !== undefined ? options.primerConc : 200e-9;
        const evalTempK = options.evalTempK !== undefined ? options.evalTempK : 310.15;
        const useOwczarzy = options.useOwczarzy !== undefined ? options.useOwczarzy : true;

        let dH = 0;
        let dS = 0;
        let gcCount = 0;

        for (let i = 0; i < seq.length - 1; i++) {
            const pair = seq.substring(i, i + 2);
            if (pair[0] === 'G' || pair[0] === 'C') gcCount++;
            if (DB.NN_WC[pair]) {
                dH += DB.NN_WC[pair].dH;
                dS += DB.NN_WC[pair].dS;
            }
        }
        if (seq[seq.length - 1] === 'G' || seq[seq.length - 1] === 'C') gcCount++;

        const fGC = gcCount / seq.length;

        dH += DB.INITIATION[seq[0]].dH + DB.INITIATION[seq[seq.length - 1]].dH;
        dS += DB.INITIATION[seq[0]].dS + DB.INITIATION[seq[seq.length - 1]].dS;

        const symmetric = isPalindrome(seq);
        if (symmetric) {
            dS -= 1.4;
        }

        const saltMod = calculateSaltCorrection(
            seq.length,
            fGC,
            naConc,
            kConc,
            trisConc,
            mgConc,
            dNTPConc,
            useOwczarzy
        );
        dS += saltMod.dS_corr;

        const effectivePrimerConc = symmetric ? primerConc : primerConc / 4;

        let tmK = 0;
        if (saltMod.tm_inv_corr !== null) {
            const baseInvTm = (dS + DB.CONSTANTS.R * Math.log(effectivePrimerConc)) / (dH * 1000);
            tmK = 1 / (baseInvTm + saltMod.tm_inv_corr);
        } else {
            tmK = (dH * 1000) / (dS + DB.CONSTANTS.R * Math.log(effectivePrimerConc));
        }

        return {
            tm: tmK - 273.15,
            dG: (dH * 1000 - evalTempK * dS) / 1000,
            dH,
            dS
        };
    }

    function getDG(dH, dS, tempK) {
        return (dH * 1000 - tempK * dS) / 1000;
    }

    function calculateGCPercent(seq) {
        if (!seq) return 0;
        const gcCount = (seq.match(/[GC]/g) || []).length;
        return (gcCount / seq.length) * 100;
    }

    function buildPrimerScoreReport(thermo, hairpin, selfDimer, metrics, seq, ta) {
        const tm = thermo?.tm ?? 0;
        const gc = calculateGCPercent(seq);
        const hairpinTm = hairpin?.Tm ?? 0;
        const dimerAnyDG = selfDimer?.ANY?.dG ?? 0;
        const dimerEndDG = selfDimer?.END?.dG ?? 0;
        const end3DG = metrics?.end3StabilityDG ?? 0;
        const gcClamp = metrics?.gcClampCount ?? 0;
        const hasRuns = !!metrics?.hasRuns;
        const seqLength = seq?.length ?? 0;

        let score = 100;
        const penalties = [];

        function addPenalty(points, label, valueText, severity) {
            score -= points;
            penalties.push({
                points,
                label,
                valueText,
                severity
            });
        }

        if (tm < 52 || tm > 68) addPenalty(15, 'Tm', `${tm.toFixed(1)} °C`, 'outside target range');
        else if (tm < 57 || tm > 63) addPenalty(5, 'Tm', `${tm.toFixed(1)} °C`, 'marginal');

        if (gc < 30 || gc > 70) addPenalty(15, 'GC%', `${gc.toFixed(1)}%`, 'outside target range');
        else if (gc < 40 || gc > 60) addPenalty(5, 'GC%', `${gc.toFixed(1)}%`, 'marginal');

        if (hairpinTm > 0 && hairpinTm > ta - 5) {
            addPenalty(20, 'Hairpin Tm', `${hairpinTm.toFixed(1)} °C (Ta = ${ta.toFixed(0)} °C)`, 'stable at annealing temp');
        } else if (hairpinTm > 0 && hairpinTm > ta - 15) {
            addPenalty(10, 'Hairpin Tm', `${hairpinTm.toFixed(1)} °C (Ta = ${ta.toFixed(0)} °C)`, 'partially stable near annealing');
        }

        if (hairpin?.structureFound && hairpin?.isThreePrimeExtensible &&
            hairpinTm > 0 && hairpinTm > ta - 15) {
            addPenalty(15, "Hairpin 3' self-priming",
                "paired 3' end with contiguous stem and 5' template overhang",
                'can self-extend');
        }

        const hasExtensibleDimer = !!selfDimer?.END?.isThreePrimeExtensible;
        const sameStructure = hasExtensibleDimer && selfDimer?.END?.sameAsGlobal === true;

        if (sameStructure) {
            // Same dimer -- only penalize as 3' extensible (the more severe interpretation)
            if (dimerEndDG < -8.0) addPenalty(25, "Self-dimer (3' extensible)", `${dimerEndDG.toFixed(1)} kcal/mol`, 'strong');
            else if (dimerEndDG < -5.0) addPenalty(15, "Self-dimer (3' extensible)", `${dimerEndDG.toFixed(1)} kcal/mol`, 'moderate');
        } else {
            // Different structures -- penalize each independently
            if (dimerAnyDG < -8.0) addPenalty(20, 'Self-dimer global ΔG', `${dimerAnyDG.toFixed(1)} kcal/mol`, 'strong');
            else if (dimerAnyDG < -5.0) addPenalty(10, 'Self-dimer global ΔG', `${dimerAnyDG.toFixed(1)} kcal/mol`, 'moderate');

            if (hasExtensibleDimer && dimerEndDG < -8.0) addPenalty(25, "3' extensible dimer ΔG", `${dimerEndDG.toFixed(1)} kcal/mol`, 'strong');
            else if (hasExtensibleDimer && dimerEndDG < -5.0) addPenalty(15, "3' extensible dimer ΔG", `${dimerEndDG.toFixed(1)} kcal/mol`, 'moderate');
        }

        if (end3DG < -10.0) addPenalty(10, "3' end stability", `${end3DG.toFixed(1)} kcal/mol`, 'over-stable');
        else if (end3DG < -9.0) addPenalty(5, "3' end stability", `${end3DG.toFixed(1)} kcal/mol`, 'marginal');

        if (gcClamp === 0 || gcClamp === 5) addPenalty(10, 'GC clamp', `${gcClamp}/5 terminal G/C`, 'poor');
        else if (gcClamp === 4) addPenalty(5, 'GC clamp', `${gcClamp}/5 terminal G/C`, 'heavy');

        if (hasRuns) addPenalty(10, 'Runs', 'detected', 'homopolymer run');

        if (seqLength < 15 || seqLength > 30) addPenalty(10, 'Length', `${seqLength} nt`, 'outside target range');
        else if (seqLength < 18 || seqLength > 25) addPenalty(5, 'Length', `${seqLength} nt`, 'marginal');

        return {
            score: Math.max(0, Math.round(score)),
            penalties
        };
    }

    function calculatePrimerScore(thermo, hairpin, selfDimer, metrics, seq, ta) {
        return buildPrimerScoreReport(thermo, hairpin, selfDimer, metrics, seq, ta).score;
    }

    function getPrimerScoreDetails(thermo, hairpin, selfDimer, metrics, seq, ta) {
        return buildPrimerScoreReport(thermo, hairpin, selfDimer, metrics, seq, ta);
    }

    function emptyHairpin() {
        return {
            dG: 0,
            dH: 0,
            dS: 0,
            Tm: 0,
            structure: null,
            structureFound: false,
            basePairCount: 0,
            loopLength: 0,
            loopSequence: '',
            specialTriloop: false,
            threePrimeRunLength: 0,
            fivePrimeTemplateOverhang: 0,
            isThreePrimeExtensible: false
        };
    }

    function emptyDimer() {
        return {
            dG: 0,
            dH: 0,
            dS: 0,
            Tm: 0,
            alignEnd1: -1,
            alignEnd2: -1,
            structure: null,
            structureFound: false,
            message: '',
            basePairCount: 0,
            wobblePairCount: 0,
            mismatchCount: 0,
            topThreePrimeRunLength: 0,
            bottomThreePrimeRunLength: 0,
            topThreePrimeExtensible: false,
            bottomThreePrimeExtensible: false,
            isThreePrimeExtensible: false,
            extensionPrimer: null,
            sameAsGlobal: false
        };
    }

    function extractPrimer3Lines(structure) {
        return structure.replace(/\n+$/, '').split('\n').map(line => {
            const tabIndex = line.indexOf('\t');
            return {
                label: tabIndex >= 0 ? line.slice(0, tabIndex) : '',
                text: tabIndex >= 0 ? line.slice(tabIndex + 1) : line
            };
        });
    }

    function overlayRows(rows) {
        const width = rows.reduce((max, row) => Math.max(max, row.length), 0);
        const overlay = new Array(width).fill(' ');

        rows.forEach(row => {
            for (let i = 0; i < row.length; i++) {
                if (row[i] !== ' ') overlay[i] = row[i];
            }
        });

        return overlay.join('');
    }

    function isBase(char) {
        return char === 'A' || char === 'T' || char === 'C' || char === 'G';
    }

    function isWatsonCrickPair(topBase, bottomBase) {
        return isBase(topBase) && isBase(bottomBase) && comp(topBase) === bottomBase;
    }

    function isWobblePair(topBase, bottomBase) {
        const pair = `${topBase}${bottomBase}`;
        return pair === 'GT' || pair === 'TG' || pair === 'GA' || pair === 'AG';
    }

    function getPairType(topBase, bottomBase) {
        if (isWatsonCrickPair(topBase, bottomBase)) return 'wc';
        if (isWobblePair(topBase, bottomBase)) return 'wobble';
        return 'mismatch';
    }

    function mapBaseIndices(chars) {
        const indices = new Array(chars.length).fill(-1);
        let baseIndex = 0;
        for (let i = 0; i < chars.length; i++) {
            if (isBase(chars[i])) {
                indices[i] = baseIndex;
                baseIndex++;
            }
        }
        return indices;
    }

    function parseDimerPairs(structure) {
        if (!structure) return [];

        const lines = extractPrimer3Lines(structure);
        const topRows = lines.filter(line => line.label === 'SEQ').map(line => line.text);
        const bottomRows = lines.filter(line => line.label === 'STR').map(line => line.text);
        if (topRows.length === 0 || bottomRows.length === 0) return [];

        const topRaw = overlayRows(topRows).replace(/-/g, ' ');
        const bottomRaw = overlayRows(bottomRows).replace(/-/g, ' ');
        const width = Math.max(topRaw.length, bottomRaw.length);
        const topChars = topRaw.padEnd(width, ' ').split('');
        const bottomChars = bottomRaw.padEnd(width, ' ').split('');
        const topIndices = mapBaseIndices(topChars);
        const bottomIndices = mapBaseIndices(bottomChars);
        const pairs = [];

        for (let i = 0; i < width; i++) {
            if (!isBase(topChars[i]) || !isBase(bottomChars[i])) continue;
            pairs.push({
                col: i,
                topIndex: topIndices[i],
                bottomIndex: bottomIndices[i],
                topBase: topChars[i],
                bottomBase: bottomChars[i],
                type: getPairType(topChars[i], bottomChars[i])
            });
        }

        return pairs;
    }

    const MIN_THREE_PRIME_TEMPLATE_OVERHANG = 2;
    const MIN_THREE_PRIME_RUN_PAIRS = 2;

    function pairsBelongToSameRun(leftPair, rightPair) {
        return rightPair.col - leftPair.col === 1 &&
               rightPair.topIndex - leftPair.topIndex === 1 &&
               rightPair.bottomIndex - leftPair.bottomIndex === 1 &&
               leftPair.type === 'wc' && rightPair.type === 'wc';
    }

    function getTopThreePrimeExtensionInfo(pairs, seq1Length, seq2Length) {
        const canonicalPairs = pairs.filter(pair => pair.type === 'wc');
        if (canonicalPairs.length === 0) return null;

        const terminalIndex = seq1Length - 1;
        const terminalPairIndex = canonicalPairs.findIndex(pair => pair.topIndex === terminalIndex);
        if (terminalPairIndex === -1) return null;

        let runStart = terminalPairIndex;
        while (runStart > 0 && pairsBelongToSameRun(canonicalPairs[runStart - 1], canonicalPairs[runStart])) {
            runStart--;
        }

        const terminalPair = canonicalPairs[terminalPairIndex];
        return {
            runLength: terminalPairIndex - runStart + 1,
            templateOverhang: (seq2Length - 1) - terminalPair.bottomIndex
        };
    }

    function getBottomThreePrimeExtensionInfo(pairs) {
        const canonicalPairs = pairs.filter(pair => pair.type === 'wc');
        if (canonicalPairs.length === 0) return null;

        const terminalPairIndex = canonicalPairs.findIndex(pair => pair.bottomIndex === 0);
        if (terminalPairIndex === -1) return null;

        let runEnd = terminalPairIndex;
        while (runEnd < canonicalPairs.length - 1 &&
               pairsBelongToSameRun(canonicalPairs[runEnd], canonicalPairs[runEnd + 1])) {
            runEnd++;
        }

        const terminalPair = canonicalPairs[terminalPairIndex];
        return {
            runLength: runEnd - terminalPairIndex + 1,
            templateOverhang: terminalPair.topIndex
        };
    }

    function analyzeDimerStructure(structure, seq1Length, seq2Length) {
        const pairs = parseDimerPairs(structure);
        const topExtension = getTopThreePrimeExtensionInfo(pairs, seq1Length, seq2Length);
        const bottomExtension = getBottomThreePrimeExtensionInfo(pairs);
        const topThreePrimeExtensible = !!topExtension &&
            topExtension.runLength >= MIN_THREE_PRIME_RUN_PAIRS &&
            topExtension.templateOverhang >= MIN_THREE_PRIME_TEMPLATE_OVERHANG;
        const bottomThreePrimeExtensible = !!bottomExtension &&
            bottomExtension.runLength >= MIN_THREE_PRIME_RUN_PAIRS &&
            bottomExtension.templateOverhang >= MIN_THREE_PRIME_TEMPLATE_OVERHANG;

        return {
            basePairCount: pairs.filter(pair => pair.type === 'wc').length,
            wobblePairCount: pairs.filter(pair => pair.type === 'wobble').length,
            mismatchCount: pairs.filter(pair => pair.type === 'mismatch').length,
            topThreePrimeRunLength: topExtension?.runLength ?? 0,
            bottomThreePrimeRunLength: bottomExtension?.runLength ?? 0,
            topTemplateOverhang: topExtension?.templateOverhang ?? 0,
            bottomTemplateOverhang: bottomExtension?.templateOverhang ?? 0,
            topThreePrimeExtensible,
            bottomThreePrimeExtensible,
            isThreePrimeExtensible: topThreePrimeExtensible || bottomThreePrimeExtensible
        };
    }

    function analyzeHairpinStructure(structure) {
        const empty = {
            basePairCount: 0,
            loopLength: 0,
            loopSequence: '',
            specialTriloop: false,
            threePrimeRunLength: 0,
            fivePrimeTemplateOverhang: 0,
            isThreePrimeExtensible: false
        };
        if (!structure) return empty;

        const lines = extractPrimer3Lines(structure);
        const foldLine = lines.find(line => line.label === 'SEQ')?.text || '';
        const seqLine = lines.find(line => line.label === 'STR')?.text || '';
        if (!foldLine || !seqLine) return empty;

        const leftPairs = [];
        const rightPairs = [];
        for (let i = 0; i < foldLine.length; i++) {
            if (foldLine[i] === '/') leftPairs.push(i);
            if (foldLine[i] === '\\') rightPairs.push(i);
        }

        const basePairCount = Math.min(leftPairs.length, rightPairs.length);
        if (basePairCount === 0) return empty;

        const innerLeft = leftPairs[leftPairs.length - 1];
        const innerRight = rightPairs[0];
        const loopSequence = seqLine
            .slice(innerLeft + 1, innerRight)
            .split('')
            .filter(isBase)
            .join('');
        const specialTriloop = loopSequence.length === 3 &&
            loopSequence[0] === 'G' && loopSequence[2] === 'A' &&
            isWatsonCrickPair(seqLine[innerLeft], seqLine[innerRight]);

        let lastBaseColumn = -1;
        for (let i = seqLine.length - 1; i >= 0; i--) {
            if (isBase(seqLine[i])) {
                lastBaseColumn = i;
                break;
            }
        }

        const outerLeft = leftPairs[0];
        const outerRight = rightPairs[rightPairs.length - 1];
        const threePrimeEndPaired = lastBaseColumn === outerRight;
        let threePrimeRunLength = 0;

        if (threePrimeEndPaired) {
            for (let offset = 0; offset < basePairCount; offset++) {
                const left = leftPairs[offset];
                const right = rightPairs[rightPairs.length - 1 - offset];
                if (offset > 0) {
                    const previousLeft = leftPairs[offset - 1];
                    const previousRight = rightPairs[rightPairs.length - offset];
                    if (left - previousLeft !== 1 || previousRight - right !== 1) break;
                }
                if (!isWatsonCrickPair(seqLine[left], seqLine[right])) break;
                threePrimeRunLength++;
            }
        }

        const fivePrimeTemplateOverhang = seqLine
            .slice(0, outerLeft)
            .split('')
            .filter(isBase)
            .length;
        const isThreePrimeExtensible = threePrimeRunLength >= MIN_THREE_PRIME_RUN_PAIRS &&
            fivePrimeTemplateOverhang >= MIN_THREE_PRIME_TEMPLATE_OVERHANG;

        return {
            basePairCount,
            loopLength: loopSequence.length,
            loopSequence,
            specialTriloop,
            threePrimeRunLength,
            fivePrimeTemplateOverhang,
            isThreePrimeExtensible
        };
    }

    function normalizeDimerResult(result) {
        return result?.structureFound && result.structure ? result : emptyDimer();
    }

    function annotateDimerResult(result, seq1Length, seq2Length) {
        const normalized = normalizeDimerResult(result);
        if (!normalized.structureFound) return normalized;
        return {
            ...normalized,
            ...analyzeDimerStructure(normalized.structure, seq1Length, seq2Length)
        };
    }

    function normalizedStructure(structure) {
        return (structure || '')
            .replace(/[ \t]+$/gm, '')
            .replace(/\n+$/, '');
    }

    function chooseMostStable(candidates) {
        return candidates.reduce((best, candidate) => (
            !best || candidate.dG < best.dG ? candidate : best
        ), null);
    }

    function primer3Ready() {
        return !!(window.Primer3 && typeof window.Primer3.calcHairpin === 'function' && window.Primer3.isReady());
    }

    function findBestHairpin(seq, options = {}) {
        seq = seq.toUpperCase();
        if (seq.length < 5 || !primer3Ready()) return emptyHairpin();

        const result = window.Primer3.calcHairpin(seq, options);
        if (!result.structureFound) return emptyHairpin();
        return {
            ...result,
            ...analyzeHairpinStructure(result.structure)
        };
    }

    function findBestDimer(p1, p2, options = {}) {
        p1 = p1.toUpperCase();
        p2 = p2.toUpperCase();

        if (!p1 || !p2 || !primer3Ready()) {
            return { ANY: emptyDimer(), END: emptyDimer() };
        }

        const any = annotateDimerResult(
            window.Primer3.calcDimer(p1, p2, options, 'ANY'),
            p1.length,
            p2.length
        );
        const endCandidates = [];
        const end1 = annotateDimerResult(
            window.Primer3.calcDimer(p1, p2, options, 'END'),
            p1.length,
            p2.length
        );

        if (end1.structureFound &&
            (end1.topThreePrimeExtensible || (p1 === p2 && end1.bottomThreePrimeExtensible))) {
            endCandidates.push({ ...end1, extensionPrimer: 1 });
        }

        if (p1 !== p2) {
            const end2 = annotateDimerResult(
                window.Primer3.calcDimer(p2, p1, options, 'END'),
                p2.length,
                p1.length
            );
            if (end2.structureFound && end2.topThreePrimeExtensible) {
                endCandidates.push({ ...end2, extensionPrimer: 2 });
            }
        }

        if (any.structureFound && any.topThreePrimeExtensible) {
            endCandidates.push({ ...any, extensionPrimer: 1, sameAsGlobal: true });
        }
        if (any.structureFound && any.bottomThreePrimeExtensible) {
            endCandidates.push({
                ...any,
                extensionPrimer: p1 === p2 ? 1 : 2,
                sameAsGlobal: true
            });
        }

        let end = chooseMostStable(endCandidates);
        if (!end) {
            end = {
                ...emptyDimer(),
                message: "No 3' extensible dimer found"
            };
        } else {
            const sameAsGlobal = end.sameAsGlobal ||
                normalizedStructure(end.structure) === normalizedStructure(any.structure);
            end = {
                ...end,
                isThreePrimeExtensible: true,
                sameAsGlobal
            };
        }

        return {
            ANY: any,
            END: end
        };
    }

    window.BioAlgorithms = {
        calculateSaltCorrection,
        calculateThermodynamics,
        calculatePrimerScore,
        findBestDimer,
        findBestHairpin,
        analyzeDimerStructure,
        analyzeHairpinStructure,
        getDG,
        getPrimerScoreDetails,
        comp,
        isPalindrome
    };
})();
