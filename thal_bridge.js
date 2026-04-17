(function() {
    const root = typeof window !== 'undefined' ? window : globalThis;
    let modulePromise = null;
    let thalModule = null;
    let thalReady = false;
    let initFailed = false;

    function wasmPath(path) {
        return `primer3_wasm/dist/${path}`;
    }

    function readLocalWasmBinary(path) {
        if (typeof readFile !== 'function') return null;

        let raw;
        try {
            raw = readFile(path);
        } catch (_) {
            return null;
        }
        if (raw instanceof Uint8Array) return raw;
        if (raw instanceof ArrayBuffer) return new Uint8Array(raw);
        if (typeof raw !== 'string') return null;

        const bytes = new Uint8Array(raw.length);
        for (let i = 0; i < raw.length; i++) {
            bytes[i] = raw.charCodeAt(i) & 0xff;
        }
        return bytes;
    }

    function buildModuleConfig() {
        const config = {
            locateFile(path) {
                return wasmPath(path);
            }
        };

        const wasmBinary = readLocalWasmBinary(wasmPath('thal.wasm'));
        if (wasmBinary) {
            config.wasmBinary = wasmBinary;
        }

        return config;
    }

    function optionsToThalArgs(options = {}) {
        return {
            mv: ((options.naConc ?? 0.05) + (options.kConc ?? 0) + ((options.trisConc ?? 0) / 2)) * 1000,
            dv: (options.mgConc ?? 0.0015) * 1000,
            dntp: (options.dNTPConc ?? 0.0008) * 1000,
            dna: (options.primerConc ?? 250e-9) * 1e9,
            tempC: (options.evalTempK ?? 310.15) - 273.15,
            maxLoop: options.maxLoopSize ?? 30
        };
    }

    async function initThal() {
        if (thalReady) return true;
        if (initFailed) return false;
        if (!modulePromise) {
            const loader = root.ThalModule || globalThis.ThalModule;
            if (typeof loader !== 'function') {
                console.error('Primer3 thal loader is missing.');
                initFailed = true;
                return false;
            }
            modulePromise = loader(buildModuleConfig());
        }

        thalModule = await modulePromise;
        const initResult = thalModule.ccall('thal_init', 'number', [], []);
        if (initResult !== 0) {
            console.error('Primer3 thal initialization failed:', initResult);
            initFailed = true;
            return false;
        }

        thalReady = true;
        return true;
    }

    function calcHairpin(seq, options = {}) {
        if (!thalReady) return { dG: 0, Tm: 0, structure: null, structureFound: false };

        const args = optionsToThalArgs(options);
        thalModule.ccall(
            'thal_calc_hairpin',
            'number',
            ['string', 'number', 'number', 'number', 'number', 'number', 'number'],
            [seq, args.mv, args.dv, args.dntp, args.dna, args.tempC, args.maxLoop]
        );

        const found = thalModule.ccall('thal_get_structure_found', 'number', [], []) === 1;
        return {
            dG: thalModule.ccall('thal_get_dG', 'number', [], []) / 1000,
            dH: thalModule.ccall('thal_get_dH', 'number', [], []),
            dS: thalModule.ccall('thal_get_dS', 'number', [], []),
            Tm: thalModule.ccall('thal_get_tm', 'number', [], []),
            structure: found ? thalModule.UTF8ToString(thalModule.ccall('thal_get_structure', 'number', [], [])) : null,
            structureFound: found
        };
    }

    function calcDimer(seq1, seq2, options = {}, mode = 'ANY') {
        if (!thalReady) return { dG: 0, Tm: 0, structure: null, structureFound: false };

        const args = optionsToThalArgs(options);
        const modeInt = mode === 'END' ? 2 : 1;
        thalModule.ccall(
            'thal_calc_dimer',
            'number',
            ['string', 'string', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
            [seq1, seq2, args.mv, args.dv, args.dntp, args.dna, args.tempC, args.maxLoop, modeInt]
        );

        const found = thalModule.ccall('thal_get_structure_found', 'number', [], []) === 1;
        return {
            dG: thalModule.ccall('thal_get_dG', 'number', [], []) / 1000,
            dH: thalModule.ccall('thal_get_dH', 'number', [], []),
            dS: thalModule.ccall('thal_get_dS', 'number', [], []),
            Tm: thalModule.ccall('thal_get_tm', 'number', [], []),
            alignEnd1: thalModule.ccall('thal_get_align_end_1', 'number', [], []),
            alignEnd2: thalModule.ccall('thal_get_align_end_2', 'number', [], []),
            structure: found ? thalModule.UTF8ToString(thalModule.ccall('thal_get_structure', 'number', [], [])) : null,
            structureFound: found
        };
    }

    root.Primer3 = {
        init: initThal,
        calcHairpin,
        calcDimer,
        isReady: () => thalReady
    };
})();
