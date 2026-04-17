document.addEventListener('DOMContentLoaded', async () => {
    const primer1Input = document.getElementById('primer1');
    const primer2Input = document.getElementById('primer2');
    const clearPrimersBtn = document.getElementById('clearPrimersBtn');
    const sequenceLayout = document.querySelector('.sequence-layout');
    const analyzeBtn = document.getElementById('analyzeBtn');
    const resultsContainer = document.getElementById('resultsContainer');

    analyzeBtn.disabled = true;
    analyzeBtn.textContent = 'Loading engine...';

    const primer3Ready = window.Primer3 && await window.Primer3.init();
    if (!primer3Ready) {
        analyzeBtn.textContent = 'Engine failed to load';
        return;
    }

    analyzeBtn.disabled = false;
    analyzeBtn.textContent = 'Analyze Primers';

    if (clearPrimersBtn && sequenceLayout) {
        clearPrimersBtn.addEventListener('click', () => {
            primer1Input.value = '';
            primer2Input.value = '';
            primer1Input.focus();
            syncClearButtonGeometry();
        });

        window.addEventListener('resize', syncClearButtonGeometry);
        requestAnimationFrame(syncClearButtonGeometry);

        if (document.fonts?.ready) {
            document.fonts.ready.then(syncClearButtonGeometry).catch(() => {});
        }
    }

    analyzeBtn.addEventListener('click', () => {
        const p1 = primer1Input.value.replace(/\s+/g, '').toUpperCase();
        const p2 = primer2Input.value.replace(/\s+/g, '').toUpperCase();

        if (!isValidDNA(p1)) {
            alert('Primer 1 contains invalid characters.');
            return;
        }
        if (p2 && !isValidDNA(p2)) {
            alert('Primer 2 contains invalid characters.');
            return;
        }
        if (!p1) {
            alert('Please enter at least Primer 1.');
            return;
        }

        const options = getOptions();
        const structureOptions = {
            ...options,
            evalTempK: options.ta + 273.15
        };
        const thresholds = getThresholds();

        if (options.dNTPConc >= options.mgConc) {
            alert("WARNING: dNTP concentration is >= Mg2+ concentration. Free Mg2+ is effectively zero. Tm predictions may be unreliable.");
        }

        resultsContainer.classList.remove('hidden');
        document.getElementById('p1Results').classList.remove('hidden');
        document.getElementById('secondaryStructuresGrid').classList.remove('hidden');

        const hairpinContainer = document.getElementById('hairpinResults');
        hairpinContainer.innerHTML = '<h2 class="primer-title" style="color: var(--accent-rose);">Hairpins</h2>';
        hairpinContainer.classList.remove('hidden');

        const p1Data = renderPrimer('p1Results', 'Primer 1', p1, options, structureOptions, thresholds);
        appendHairpinToContainer(hairpinContainer, 'Primer 1 Hairpin', p1Data.hairpin, thresholds, options.ta);

        if (p2) {
            document.getElementById('p2Results').classList.remove('hidden');
            document.getElementById('crossResults').classList.remove('hidden');

            const p2Data = renderPrimer('p2Results', 'Primer 2', p2, options, structureOptions, thresholds);
            appendHairpinToContainer(hairpinContainer, 'Primer 2 Hairpin', p2Data.hairpin, thresholds, options.ta);

            const cross = window.BioAlgorithms.findBestDimer(p1, p2, structureOptions);
            renderCrossDimer('crossResults', cross, p1Data.thermo.tm, p2Data.thermo.tm, thresholds, options.ta);
        } else {
            document.getElementById('p2Results').classList.add('hidden');
            document.getElementById('crossResults').classList.add('hidden');
        }
    });

    function syncClearButtonGeometry() {
        if (!clearPrimersBtn || !sequenceLayout) return;

        const layoutRect = sequenceLayout.getBoundingClientRect();
        const primer1Rect = primer1Input.getBoundingClientRect();
        const primer2Rect = primer2Input.getBoundingClientRect();

        const topOffset = Math.max(0, primer1Rect.top - layoutRect.top);
        const totalHeight = Math.max(0, primer2Rect.bottom - primer1Rect.top);

        clearPrimersBtn.style.marginTop = `${topOffset}px`;
        clearPrimersBtn.style.height = `${totalHeight}px`;
    }

    function isValidDNA(seq) {
        return /^[ATCG]+$/.test(seq);
    }

    function calculateGC(seq) {
        if (!seq) return 0;
        const gc = (seq.match(/[GC]/g) || []).length;
        return ((gc / seq.length) * 100).toFixed(1);
    }

    function getOptions() {
        return {
            naConc: parseFloat(document.getElementById('p_na').value) / 1000,
            mgConc: parseFloat(document.getElementById('p_mg').value) / 1000,
            dNTPConc: parseFloat(document.getElementById('p_dntp').value) / 1000,
            kConc: parseFloat(document.getElementById('p_k').value) / 1000,
            trisConc: parseFloat(document.getElementById('p_tris').value) / 1000,
            primerConc: parseFloat(document.getElementById('p_primer_conc').value) * 1e-9,
            evalTempK: parseFloat(document.getElementById('p_eval_t').value) + 273.15,
            ta: parseFloat(document.getElementById('p_ta').value),
            useOwczarzy: document.getElementById('p_tm_formula').value === 'owczarzy',
            maxLoopSize: parseInt(document.getElementById('p_max_hairpin').value, 10)
        };
    }

    function getThresholds() {
        return {
            maxHairpinDG: parseFloat(document.getElementById('t_hairpin_dg').value),
            maxDimerAnyDG: parseFloat(document.getElementById('t_dimer_any').value),
            maxDimerEndDG: parseFloat(document.getElementById('t_dimer_end').value),
            maxTmDiff: parseFloat(document.getElementById('t_tm_diff').value),
            maxEnd3StabDG: parseFloat(document.getElementById('t_end3_stab').value)
        };
    }

    function analyzeSequenceMetrics(seq, options) {
        seq = seq.toUpperCase();
        const metrics = {
            gcClampCount: 0,
            hasRuns: false,
            hasRepeats: false,
            end3StabilityDG: 0,
            end5StabilityDG: 0
        };
        if (seq.length < 5) return metrics;

        const last5 = seq.substring(seq.length - 5);
        const first5 = seq.substring(0, 5);
        metrics.gcClampCount = (last5.match(/[GC]/g) || []).length;
        metrics.end3StabilityDG = window.BioAlgorithms.calculateThermodynamics(last5, options).dG;
        metrics.end5StabilityDG = window.BioAlgorithms.calculateThermodynamics(first5, options).dG;
        metrics.hasRuns = /(.)\1{3,}/.test(seq);
        metrics.hasRepeats = /(..)\1{3,}/.test(seq);
        return metrics;
    }

    function renderPrimer(containerId, title, seq, thermoOptions, structureOptions, thresholds) {
        const container = document.getElementById(containerId);
        container.innerHTML = '';

        const template = document.getElementById('primerReportTemplate');
        const content = template.content.cloneNode(true);

        content.querySelector('.primer-title').textContent = title;
        content.querySelector('.sequence-display').textContent = `5'-${seq}-3'`;

        const thermo = window.BioAlgorithms.calculateThermodynamics(seq, thermoOptions);
        const templateHybrid = window.BioAlgorithms.calculateThermodynamics(seq, structureOptions);
        const selfDimer = window.BioAlgorithms.findBestDimer(seq, seq, structureOptions);
        const hairpin = window.BioAlgorithms.findBestHairpin(seq, structureOptions);
        const meta = analyzeSequenceMetrics(seq, thermoOptions);

        const tmCell = content.querySelector('.t-tm');
        tmCell.textContent = thermo.tm > 0 ? `${thermo.tm.toFixed(1)} °C` : '--';
        tmCell.className = `metric-value t-tm ${getTrafficColor(thermo.tm, 57, 63, 52, 68, true)}`;

        const templateDGCell = content.querySelector('.t-template-dg');
        templateDGCell.textContent = `${templateHybrid.dG.toFixed(2)} kcal/mol`;
        templateDGCell.className = 'metric-value t-template-dg status-default';

        const gcCell = content.querySelector('.t-gc');
        const gc = parseFloat(calculateGC(seq));
        gcCell.textContent = `${gc}%`;
        gcCell.className = `metric-value t-gc ${getTrafficColor(gc, 40, 60, 30, 70, true)}`;

        content.querySelector('.t-len').textContent = `${seq.length} nt`;

        const score = window.BioAlgorithms.calculatePrimerScore(thermo, hairpin, selfDimer, meta, seq, thermoOptions.ta);
        const scoreReport = window.BioAlgorithms.getPrimerScoreDetails(thermo, hairpin, selfDimer, meta, seq, thermoOptions.ta);
        const scoreCell = content.querySelector('.t-score');
        scoreCell.textContent = `${score}`;
        scoreCell.className = `metric-value t-score ${scoreAlert(score)}`;
        populateScoreBreakdown(content.querySelector('.metric-breakdown'), scoreReport.penalties);

        const structContainer = content.querySelector('.structures-container');
        const structureTempC = structureOptions.evalTempK - 273.15;

        const sdAnyAlert = structureAlert(selfDimer.ANY.dG, thresholds.maxDimerAnyDG);
        if (selfDimer.ANY.structureFound) {
            structContainer.appendChild(createAlignmentBlock('Self-Dimer (Global)', selfDimer.ANY, sdAnyAlert, structureTempC, true));
        } else {
            structContainer.appendChild(createStructureLine('Self-Dimer (Global)', selfDimer.ANY.dG, structureTempC, sdAnyAlert, selfDimer.ANY.message));
        }

        const sdEndAlert = structureAlert(selfDimer.END.dG, thresholds.maxDimerEndDG);
        if (selfDimer.END.structureFound) {
            structContainer.appendChild(createAlignmentBlock("Self-Dimer (3' Extensible)", selfDimer.END, sdEndAlert, structureTempC, true));
        } else {
            structContainer.appendChild(createStructureLine("Self-Dimer (3' Extensible)", selfDimer.END.dG, structureTempC, sdEndAlert, selfDimer.END.message));
        }

        const checks = content.querySelector('.checks-container');
        const end3Class = meta.end3StabilityDG < thresholds.maxEnd3StabDG ? 'status-red' : 'status-green';
        checks.appendChild(createCheckLine("3' End Stability (last 5)", `${meta.end3StabilityDG.toFixed(2)} kcal/mol`, end3Class));

        const gcClass = meta.gcClampCount >= 1 && meta.gcClampCount <= 3 ? 'status-green' : 'status-yellow';
        checks.appendChild(createCheckLine('GC Clamp (last 5)', `${meta.gcClampCount}/5 terminal G/C`, gcClass));

        checks.appendChild(createCheckLine('Runs (>=4)', meta.hasRuns ? 'DETECTED' : 'None', meta.hasRuns ? 'status-red' : 'status-green'));
        checks.appendChild(createCheckLine('Repeats', meta.hasRepeats ? 'DETECTED' : 'None', meta.hasRepeats ? 'status-red' : 'status-green'));

        container.appendChild(content);
        return { hairpin, thermo };
    }

    function appendHairpinToContainer(container, title, hairpin, thresholds, structureTempC) {
        const hpAlert = structureAlert(hairpin.dG, thresholds.maxHairpinDG);
        if (hairpin.structureFound) {
            container.appendChild(createHairpinBlock(title, hairpin, hpAlert, structureTempC));
        } else {
            container.appendChild(createStructureLine(title, hairpin.dG, structureTempC, hpAlert));
        }
    }

    function renderCrossDimer(containerId, dimer, p1Tm, p2Tm, thresholds, structureTempC) {
        const container = document.getElementById(containerId);
        container.innerHTML = '<h2 class="primer-title" style="color: var(--accent-rose);">Cross-Primer Dimer</h2>';
        const structContainer = document.createElement('div');

        const tmDiff = Math.abs(p1Tm - p2Tm);
        structContainer.appendChild(
            createCheckLine('Primer Pair ΔTm', `${tmDiff.toFixed(1)} °C`, tmDiffAlert(tmDiff, thresholds.maxTmDiff))
        );

        const crossAnyAlert = structureAlert(dimer.ANY.dG, thresholds.maxDimerAnyDG);
        if (dimer.ANY.structureFound) {
            structContainer.appendChild(createAlignmentBlock('Cross-Dimer (Global)', dimer.ANY, crossAnyAlert, structureTempC, true));
        } else {
            structContainer.appendChild(createStructureLine('Cross-Dimer (Global)', dimer.ANY.dG, structureTempC, crossAnyAlert, dimer.ANY.message));
        }

        const crossEndAlert = structureAlert(dimer.END.dG, thresholds.maxDimerEndDG);
        if (dimer.END.structureFound) {
            structContainer.appendChild(createAlignmentBlock("Cross-Dimer (3' Extensible)", dimer.END, crossEndAlert, structureTempC, true));
        } else {
            structContainer.appendChild(createStructureLine("Cross-Dimer (3' Extensible)", dimer.END.dG, structureTempC, crossEndAlert, dimer.END.message));
        }

        container.appendChild(structContainer);
    }

    function structureAlert(dG, threshold) {
        if (dG < threshold) return 'status-red';
        if (dG < threshold + 1) return 'status-yellow';
        return 'status-green';
    }

    function tmDiffAlert(tmDiff, threshold) {
        if (tmDiff > threshold + 1) return 'status-red';
        if (tmDiff > threshold) return 'status-yellow';
        return 'status-green';
    }

    function scoreAlert(score) {
        if (score >= 80) return 'status-green';
        if (score >= 60) return 'status-yellow';
        return 'status-red';
    }

    function getTrafficColor(val, greenMin, greenMax, yellowMin, yellowMax) {
        if (val >= greenMin && val <= greenMax) return 'status-green';
        if (val >= yellowMin && val <= yellowMax) return 'status-yellow';
        return 'status-red';
    }

    function populateScoreBreakdown(detailsEl, penalties) {
        if (!detailsEl) return;

        const list = detailsEl.querySelector('.score-penalties');
        const summary = detailsEl.querySelector('summary');
        if (!list || !summary) return;

        list.innerHTML = '';

        if (!penalties || penalties.length === 0) {
            summary.textContent = 'Penalties';
            const item = document.createElement('li');
            item.className = 'score-penalty-empty';
            item.textContent = 'No penalties applied.';
            list.appendChild(item);
            detailsEl.open = false;
            return;
        }

        summary.textContent = penalties.length === 1 ? '1 penalty' : `${penalties.length} penalties`;
        penalties.forEach(penalty => {
            const item = document.createElement('li');
            item.className = 'score-penalty-item';
            item.innerHTML = `<span class="score-penalty-points">-${penalty.points}</span><span class="score-penalty-text">${penalty.label} = ${penalty.valueText} (${penalty.severity})</span>`;
            list.appendChild(item);
        });
        detailsEl.open = false;
    }

    function createCheckLine(label, value, colorClass) {
        const div = document.createElement('div');
        div.className = 'check-item';
        div.innerHTML = `<span>${label}</span> <strong class="${colorClass}">${value}</strong>`;
        return div;
    }

    function createStructureLine(label, dG, tempC, colorClass, message = '') {
        const div = document.createElement('div');
        div.className = 'check-item';
        const text = message || `ΔG @ ${tempC.toFixed(1)} °C = ${dG.toFixed(2)} kcal/mol`;
        div.innerHTML = `<span>${label}</span> <strong class="${colorClass}">${text}</strong>`;
        return div;
    }

    function createAlignmentBlock(title, result, colorClass, tempC, allowFragmented = false) {
        const wrapper = createStructureWrapper(title, result, colorClass, tempC);
        if (!result.structure) return wrapper;

        const parsed = parsePrimer3DimerStructure(result.structure, allowFragmented);
        if (!parsed) {
            appendStyledPlaceholder(
                wrapper,
                `5' structure detected 3'`,
                '   visualization unavailable',
                `3' see metrics above 5'`
            );
            return wrapper;
        }

        appendStyledAlignment(
            wrapper,
            [
                renderDimerMainRow(parsed.topRow, parsed.matchRow),
                colorizeMatchRow(parsed.matchRow),
                renderDimerMainRow(parsed.bottomRow, parsed.matchRow)
            ]
        );
        return wrapper;
    }

    function createHairpinBlock(title, result, colorClass, tempC) {
        const wrapper = createStructureWrapper(title, result, colorClass, tempC);
        if (!result.structure) return wrapper;

        const parsed = parsePrimer3HairpinStructure(result.structure);
        if (!parsed) {
            appendStyledPlaceholder(
                wrapper,
                `5' structure detected 3'`,
                '   hairpin layout unavailable',
                `3' see metrics above 5'`
            );
            return wrapper;
        }

        appendStyledAlignment(
            wrapper,
            [
                highlightRange(parsed.topRow, parsed.topStemStart, parsed.stemLength, 'match-wc'),
                colorizeMatchRow(parsed.midRow),
                highlightRange(parsed.bottomRow, parsed.bottomStemStart, parsed.stemLength, 'match-wc')
            ]
        );
        return wrapper;
    }

    function createStructureWrapper(title, result, colorClass, tempC) {
        const wrapper = document.createElement('div');
        wrapper.className = 'alignment-block';

        const header = document.createElement('div');
        header.className = 'alignment-header';
        header.innerHTML = `<span>${title}</span>${formatStructureMetrics(result, colorClass, tempC)}`;
        wrapper.appendChild(header);
        return wrapper;
    }

    function formatStructureMetrics(result, colorClass, annealingTempC) {
        const dgText = `<span class="dg-value ${colorClass}">ΔG: ${result.dG.toFixed(2)} kcal/mol</span>`;
        if (!Number.isFinite(result.Tm) || result.Tm <= 0) {
            return `<span class="structure-metrics">${dgText}</span>`;
        }

        const tmClass = structureTmAlert(result.Tm, annealingTempC);
        const tmText = `<span class="tm-value ${tmClass}">Tm: ${result.Tm.toFixed(1)} °C</span>`;
        return `<span class="structure-metrics">${dgText}${tmText}</span>`;
    }

    function structureTmAlert(structureTm, annealingTempC) {
        if (structureTm < annealingTempC - 15) return 'status-green';
        if (structureTm <= annealingTempC - 5) return 'status-yellow';
        return 'status-red';
    }

    function appendStyledAlignment(wrapper, lines) {
        const pre = document.createElement('pre');
        const code = document.createElement('code');
        code.className = 'alignment-code';
        code.innerHTML = lines.join('\n');
        pre.appendChild(code);
        wrapper.appendChild(pre);
    }

    function appendStyledPlaceholder(wrapper, topRow, midRow, bottomRow) {
        appendStyledAlignment(wrapper, [
            `<span class="alignment-note-text">${escapeChar(topRow)}</span>`,
            `<span class="alignment-note-text">${escapeChar(midRow)}</span>`,
            `<span class="alignment-note-text">${escapeChar(bottomRow)}</span>`
        ]);
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

    function isWobblePair(a, b) {
        const pair = a + b;
        return pair === 'GT' || pair === 'TG' || pair === 'GA' || pair === 'AG';
    }

    function normalizeDimerRow(row) {
        return row.replace(/-/g, ' ');
    }

    function classifyDimerColumn(topChar, bottomChar) {
        const hasTop = isBase(topChar);
        const hasBottom = isBase(bottomChar);
        if (hasTop && hasBottom) return 'pair';
        if (hasTop) return 'topOnly';
        if (hasBottom) return 'bottomOnly';
        return 'none';
    }

    function getPairType(topChar, bottomChar) {
        if (!isBase(topChar) || !isBase(bottomChar)) return null;
        if (topChar === window.BioAlgorithms.comp(bottomChar)) return 'wc';
        if (isWobblePair(topChar, bottomChar)) return 'wobble';
        return 'mismatch';
    }

    function buildPairPrefix(topChars, bottomChars) {
        const prefix = new Array(topChars.length + 1).fill(0);
        for (let i = 0; i < topChars.length; i++) {
            prefix[i + 1] = prefix[i] + (classifyDimerColumn(topChars[i], bottomChars[i]) === 'pair' ? 1 : 0);
        }
        return prefix;
    }

    function collectDimerSegments(topChars, bottomChars) {
        const segments = [];
        let start = 0;

        while (start < topChars.length) {
            const kind = classifyDimerColumn(topChars[start], bottomChars[start]);
            let end = start + 1;
            while (end < topChars.length && classifyDimerColumn(topChars[end], bottomChars[end]) === kind) {
                end++;
            }
            segments.push({ kind, start, end: end - 1 });
            start = end;
        }

        return segments;
    }

    function isInternalLoopSegment(start, end, pairPrefix) {
        const pairsBefore = pairPrefix[start];
        const pairsAfter = pairPrefix[pairPrefix.length - 1] - pairPrefix[end + 1];
        return pairsBefore > 0 && pairsAfter > 0;
    }

    function parsePrimer3DimerStructure(structure, allowFragmented = false) {
        const lines = extractPrimer3Lines(structure);
        const topRows = lines.filter(line => line.label === 'SEQ').map(line => line.text);
        const bottomRows = lines.filter(line => line.label === 'STR').map(line => line.text);
        if (topRows.length === 0 || bottomRows.length === 0) return null;

        const topRaw = normalizeDimerRow(overlayRows(topRows));
        const bottomRaw = normalizeDimerRow(overlayRows(bottomRows));
        const width = Math.max(topRaw.length, bottomRaw.length);
        const topChars = topRaw.padEnd(width, ' ').split('');
        const bottomChars = bottomRaw.padEnd(width, ' ').split('');
        const topDisplay = topChars.slice();
        const bottomDisplay = bottomChars.slice();
        const segments = collectDimerSegments(topChars, bottomChars);
        const pairPrefix = buildPairPrefix(topChars, bottomChars);
        const matchChars = new Array(width).fill(' ');
        let hasPairs = false;

        segments.forEach(segment => {
            const internalLoop = isInternalLoopSegment(segment.start, segment.end, pairPrefix);

            if (segment.kind === 'pair') {
                hasPairs = true;
                for (let col = segment.start; col <= segment.end; col++) {
                    matchChars[col] = getMatchChar(topChars[col], bottomChars[col]);
                }
                return;
            }

            if (!internalLoop) return;

            if (segment.kind === 'topOnly') {
                for (let col = segment.start; col <= segment.end; col++) {
                    if (isBase(topDisplay[col])) topDisplay[col] = topDisplay[col].toLowerCase();
                    bottomDisplay[col] = '─';
                }
                return;
            }

            if (segment.kind === 'bottomOnly') {
                for (let col = segment.start; col <= segment.end; col++) {
                    if (isBase(bottomDisplay[col])) bottomDisplay[col] = bottomDisplay[col].toLowerCase();
                    topDisplay[col] = '─';
                }
            }
        });

        if (!hasPairs) return null;
        return {
            topRow: `5' ${trimRight(topDisplay.join(''))} 3'`,
            matchRow: `   ${trimRight(matchChars.join(''))}   `,
            bottomRow: `3' ${trimRight(bottomDisplay.join(''))} 5'`
        };
    }

    function getMatchChar(topChar, bottomChar) {
        if (!isBase(topChar) || !isBase(bottomChar)) return ' ';
        if (topChar === window.BioAlgorithms.comp(bottomChar)) return '|';
        if (isWobblePair(topChar, bottomChar)) return ':';
        return '.';
    }

    function trimRight(str) {
        return str.replace(/\s+$/, '');
    }

    function parsePrimer3HairpinStructure(structure) {
        const lines = extractPrimer3Lines(structure);
        const foldLine = lines.find(line => line.label === 'SEQ')?.text || '';
        const seqLine = lines.find(line => line.label === 'STR')?.text || '';
        if (!foldLine || !seqLine) return null;

        const leftStart = foldLine.indexOf('/');
        const leftEnd = foldLine.lastIndexOf('/');
        const rightStart = foldLine.indexOf('\\');
        const rightEnd = foldLine.lastIndexOf('\\');
        if (leftStart === -1 || leftEnd === -1 || rightStart === -1 || rightEnd === -1) return null;

        const tail5 = seqLine.slice(0, leftStart);
        const stem5 = seqLine.slice(leftStart, leftEnd + 1);
        const loopSeq = seqLine.slice(leftEnd + 1, rightStart);
        const stem3 = seqLine.slice(rightStart, rightEnd + 1);
        const tail3 = seqLine.slice(rightEnd + 1);
        const bottomStem = stem3.split('').reverse().join('');
        const bottomTail = tail3.split('').reverse().join('');
        const stemLength = Math.min(stem5.length, bottomStem.length);
        const padTop = Math.max(0, bottomTail.length - tail5.length);
        const padBot = Math.max(0, tail5.length - bottomTail.length);
        const loopTop = loopSeq.slice(0, Math.floor(loopSeq.length / 2));
        const loopMid = loopSeq.length % 2 === 1 ? loopSeq[Math.floor(loopSeq.length / 2)] : '│';
        const loopBot = loopSeq.slice(Math.ceil(loopSeq.length / 2)).split('').reverse().join('');
        return {
            topRow: `5' ${' '.repeat(padTop)}${tail5}${stem5}-${loopTop}┐`,
            midRow: `   ${' '.repeat(padTop + tail5.length)}${'|'.repeat(stemLength)} ${' '.repeat(loopTop.length)}${loopMid}`,
            bottomRow: `3' ${' '.repeat(padBot)}${bottomTail}${bottomStem}-${loopBot}┘`,
            topStemStart: 3 + padTop + tail5.length,
            bottomStemStart: 3 + padBot + bottomTail.length,
            stemLength
        };
    }

    function escapeChar(char) {
        if (char === '&') return '&amp;';
        if (char === '<') return '&lt;';
        if (char === '>') return '&gt;';
        return char;
    }

    function isLoopBase(char) {
        return char === 'a' || char === 't' || char === 'c' || char === 'g';
    }

    function renderDimerMainRow(rowStr, matchStr) {
        let result = '';
        for (let i = 0; i < rowStr.length; i++) {
            const rawChar = rowStr[i];
            const char = escapeChar(rawChar);
            if (i < 3 || i >= rowStr.length - 3) {
                result += char;
                continue;
            }

            if (rawChar === '─') {
                result += `<span class="dimer-loop-bridge">${char}</span>`;
                continue;
            }

            if (isLoopBase(rawChar)) {
                result += `<span class="dimer-loop-base">${char}</span>`;
                continue;
            }

            const matchChar = matchStr[i] || ' ';
            if (matchChar === '|') {
                result += `<span class="match-wc">${char}</span>`;
            } else if (matchChar === ':') {
                result += `<span class="match-wobble">${char}</span>`;
            } else if (matchChar === '.') {
                result += `<span class="match-mm">${char}</span>`;
            } else {
                result += char;
            }
        }
        return result;
    }

    function highlightRange(rowStr, start, length, className) {
        let result = '';
        for (let i = 0; i < rowStr.length; i++) {
            const char = escapeChar(rowStr[i]);
            if (i >= start && i < start + length) {
                result += `<span class="${className}">${char}</span>`;
            } else {
                result += char;
            }
        }
        return result;
    }

    function colorizeMatchRow(rowStr) {
        return rowStr
            .split('')
            .map(char => {
                if (char === '|') return '<span class="match-wc">|</span>';
                if (char === ':') return '<span class="match-wobble">:</span>';
                if (char === '.') return '<span class="match-mm">.</span>';
                return escapeChar(char);
            })
            .join('');
    }
});
