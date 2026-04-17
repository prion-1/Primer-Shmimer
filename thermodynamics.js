window.ThermoDB = (function() {
    const CONSTANTS = {
        R: 1.9872
    };

    const NN_WC = {
        AA: { dH: -7.9, dS: -22.2 }, TT: { dH: -7.9, dS: -22.2 },
        AT: { dH: -7.2, dS: -20.4 }, TA: { dH: -7.2, dS: -21.3 },
        CA: { dH: -8.5, dS: -22.7 }, TG: { dH: -8.5, dS: -22.7 },
        GT: { dH: -8.4, dS: -22.4 }, AC: { dH: -8.4, dS: -22.4 },
        CT: { dH: -7.8, dS: -21.0 }, AG: { dH: -7.8, dS: -21.0 },
        GA: { dH: -8.2, dS: -22.2 }, TC: { dH: -8.2, dS: -22.2 },
        CG: { dH: -10.6, dS: -27.2 }, GC: { dH: -9.8, dS: -24.4 },
        GG: { dH: -8.0, dS: -19.9 }, CC: { dH: -8.0, dS: -19.9 }
    };

    const INITIATION = {
        G: { dH: 0.1, dS: -2.8 },
        C: { dH: 0.1, dS: -2.8 },
        A: { dH: 2.3, dS: 4.1 },
        T: { dH: 2.3, dS: 4.1 }
    };

    return {
        CONSTANTS,
        INITIATION,
        NN_WC
    };
})();
