#!/usr/bin/env python3
"""Two objections to the KChIP exchange inference, worked through.

OBJECTION 1: the degradation is probably at the ER, so the short KChIP half-life is biosynthetic
triage of unassembled subunits, not post-assembly exchange.

The premise is almost certainly right about location - unassembled channel subunits are cleared by
ER quality control - but it cannot produce the observed number, for two independent reasons.

(a) WHAT THE MEASUREMENT FOLLOWS. Dorrbaum's design is dynamic SILAC run BACKWARDS from a pulse
    chase: cultures were grown 18-19 DIV in light medium, switched to heavy, and the decay of the
    PRE-EXISTING LIGHT protein was fitted over t = 0, 1, 3, 7 d. The quantity followed is therefore
    the steady-state population as it stood at the medium switch, not a cohort of newly synthesised
    protein. The bias that ER triage would create - label concentrated in young, doomed protein -
    requires a short pulse of new synthesis, which is not this design.

(b) FLUX IS NOT POOL. A pool's steady-state size is flux x lifetime, so a fast-degrading pool is a
    SMALL pool however much synthesis flows through it. For an ER pool to be the 0.82-0.96 steady-
    state fraction that a 2.02 d apparent half-life demands, the synthesis excess has to be
    enormous, and the faster the ER clearance the more enormous it gets. This script computes it.
    The conclusion inverts the objection: because ER degradation is FAST, it cannot explain an
    apparent half-life as SLOW as 2 days.

OBJECTION 2: the RNAi was assessed a week after the virus went in.

This one lands, and it corrects an over-reading. Seven days is 1.5-3.5 half-lives for these
subunits, so loss of protein at that timepoint is consistent with DPP6 being required to ASSEMBLE
new complexes while pre-existing complexes simply decay at their normal rate. It does not, on its
own, show that removing DPP6 destabilises standing complexes - which is what the two-pool model was
credited with. What discriminates is the EXTENT of the loss: assembly failure alone leaves a
predictable residue, and loss beyond that residue is the part that needs active co-degradation.

Usage:  python3 er_triage_and_rnai_timing.py [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
LN2 = np.log(2)

# neuron-enriched unless noted; Kcnd2 is mixed cultures (not quantified in neuron-enriched)
HALF = {'DPP6': 4.52, 'DPP10': 4.36, 'Kv4.3 (Kcnd3)': 3.87, 'Kv4.2 (Kcnd2, mixed)': 4.81,
        'KChIP1': 2.67, 'KChIP4': 2.02}
T_COMPLEX = 3.87          # d, Kv4.3 - the within-preparation reference
RNAI_DAYS = 7.0
KD_DELAY = [0.0, 1.0, 2.0, 3.0]   # d before knockdown is effectively complete


def main():
    # ---------------------------------------------------------------- objection 1b
    print('=' * 100)
    print('OBJECTION 1: could ER triage of unassembled KChIP produce a 2.02 d apparent half-life?')
    print('=' * 100)
    print('  A pool at steady state has size = flux x lifetime, so:')
    print('      f_free / (1 - f_free)  =  [phi / (1 - phi)] x [t_free / t_complex]')
    print('  where phi is the fraction of SYNTHESIS destroyed before assembly. Inverting it gives the')
    print('  synthesis excess each scenario demands.\n')
    print(f'  {"ER t1/2":>10} {"f_free needed":>15} {"phi (synth destroyed)":>23} '
          f'{"excess synthesis":>18}')
    rows = []
    # f_free values that a 2.02 d apparent half-life requires, from kchip_turnover_model.py
    need = {2 / 24: 0.96, 6 / 24: 0.93, 24 / 24: 0.82}
    for t_free, f in need.items():
        odds = (f / (1 - f)) * (T_COMPLEX / t_free)
        phi = odds / (1 + odds)
        label = f'{t_free*24:.0f} h'
        print(f'  {label:>10} {f:15.2f} {phi:23.4f} {1/(1-phi):17.0f}x')
        rows.append(dict(er_half_life_h=t_free * 24, f_free_needed=f, phi_synthesis_destroyed=phi,
                         fold_excess_synthesis=1 / (1 - phi)))
    print('\n  So a 2 h ER clearance would need 99.9% of all KChIP4 synthesis destroyed - roughly a')
    print('  thousandfold excess - and even a generous 24 h ER residence needs 95%. The objection')
    print('  defeats itself: ER quality control is FAST, fast pools are SMALL at steady state, and a')
    print('  small pool cannot set the decay rate of the bulk. An apparent half-life of 2.02 d is too')
    print('  SLOW to be a biosynthetic-triage signature.')
    print('\n  Independently: the measurement followed the decay of PRE-EXISTING protein after a')
    print('  light-to-heavy medium switch at 18-19 DIV, not the fate of a labelled new cohort, so the')
    print('  young-protein bias that ER triage would create is not in this design.')

    # ---------------------------------------------------------------- objection 2
    print('\n' + '=' * 100)
    print(f'OBJECTION 2: what does the {RNAI_DAYS:.0f}-day RNAi timepoint actually discriminate?')
    print('=' * 100)
    print('  Null model: DPP6 knockdown blocks assembly of NEW complexes; everything already made')
    print('  decays at its own normal rate. Protein remaining when assayed:\n')
    print(f'  {"subunit":22} ' + ' '.join(f'{"delay "+str(int(d))+"d":>11}' for d in KD_DELAY))
    res = []
    for name, th in HALF.items():
        k = LN2 / th
        cells = []
        for d in KD_DELAY:
            eff = max(RNAI_DAYS - d, 0.0)
            cells.append(np.exp(-k * eff))
        print(f'  {name:22} ' + ' '.join(f'{c*100:10.0f}%' for c in cells))
        res.append(dict(subunit=name, half_life_d=th,
                        **{f'remaining_delay{int(d)}d': round(float(np.exp(-LN2/th*max(RNAI_DAYS-d,0))), 4)
                           for d in KD_DELAY}))
    print('\n  Reading this. With no knockdown delay, assembly failure alone leaves roughly a THIRD of')
    print('  Kv4.2 and Kv4.3 at one week, and 9-16% of the KChIPs. Allowing 2-3 days for the virus to')
    print('  take hold and for existing DPP6 protein to fall, the predicted residue rises to about')
    print('  half for the Kv4 subunits.')
    print('\n  THEREFORE: the extent of the loss is the discriminator, not the fact of it.')
    print('   - If the blots showed roughly 30-50% of Kv4.2/Kv4.3 remaining, assembly failure alone')
    print('     accounts for it and the experiment says DPP6 is needed to BUILD complexes. It says')
    print('     nothing about whether standing complexes are destabilised.')
    print('   - If the blots showed near-total loss (say under 10%), that is faster than any of these')
    print('     subunits can decay on its own, and active co-degradation is required.')
    print('\n  CORRECTION TO THE RECORD: this analysis previously credited the DPP6 RNAi result with')
    print('  establishing that free subunits are degraded, i.e. with supplying the premise the')
    print('  two-pool model needs. At a single 7-day timepoint it does not establish that. Whether it')
    print('  does depends on the quantified extent of loss in the lab\'s own blots.')
    print('\n  WHAT WOULD SETTLE IT CHEAPLY: a 2-3 day timepoint. At 2 days, assembly failure predicts')
    for d in (1, 2, 3):
        k = LN2 / T_COMPLEX
        print(f'     day {d}: {np.exp(-k*d)*100:.0f}% of Kv4.3 remaining under assembly failure alone,')
    print('   so an early, steep loss of Kv4 would separate destabilisation from assembly failure')
    print('   far better than the one-week endpoint can.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    pd.DataFrame(rows).to_csv(f'{OUT}/tables/er_triage_requirement.csv', index=False)
    pd.DataFrame(res).to_csv(f'{OUT}/tables/rnai_timing_null.csv', index=False)
    print(f'\nwrote {OUT}/tables/{{er_triage_requirement,rnai_timing_null}}.csv')

    print('\n' + '=' * 100)
    print('WHERE THE EXCHANGE INFERENCE STANDS AFTER BOTH OBJECTIONS')
    print('=' * 100)
    print('  Still standing, because it never depended on the RNAi: within one preparation, KChIP4')
    print(f'  (2.02 d) and KChIP1 (2.67 d) turn over faster than Kv4.3 ({T_COMPLEX} d), and a mixture')
    print('  whose protected component decays at the slower rate cannot present the faster apparent')
    print('  rate. ER triage cannot rescue a protected-pool model, per objection 1.')
    print('  Weakened: the claim that the lab\'s RNAi independently demonstrates degradation of free')
    print('  subunits. That now rests on how complete the loss was.')
    print('  Unchanged: the co-IP test is the decisive one. What fraction of KChIP4 is bound to Kv4?')


if __name__ == '__main__':
    main()
