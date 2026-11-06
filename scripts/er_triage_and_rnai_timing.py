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

OBJECTION 2: the RNAi timing. RETRACTED SECTION - see the note below.

The RNAi design is: neurons plated, infected at DIV 2 with DPP6 RNAi or scrambled control, western
blotted at DIV 10. The comparison is knockdown vs control at ONE timepoint, with no measurement at
DIV 2. A freshly plated neuron has barely begun to express the complex, so across DIV 2-10 the pool
is being BUILT, not decaying from a plateau. The decay-null calculation this script used to perform
("protein remaining at 7 days") computed a quantity the experiment does not contain, and has been
removed rather than corrected.

What the experiment establishes instead: knocking down DPP6 while the complex is first assembled
prevents Kv4.2, Kv4.3 and KChIP3 from accumulating. Since only DPP6 is targeted, the others' newly
made protein must be failing to survive - the ER quality control of objection 1, operating during
biogenesis. So DPP6 is required for assembly and sits upstream of the rest, which agrees with its
being the longest-lived subunit in the mature complex. Neither fact bears on whether mature KChIP
exchanges, so the RNAi is out of the turnover argument entirely.

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

    print('\n' + '=' * 100)
    print('OBJECTION 2: RETRACTED - the RNAi design cannot be read as a decay experiment')
    print('=' * 100)
    print('  Infection at DIV 2, blot at DIV 10, knockdown vs scrambled control at that one')
    print('  timepoint, with no measurement at DIV 2. The complex is being BUILT across that window,')
    print('  not decaying from a plateau, so a "protein remaining" calculation computes a quantity')
    print('  the experiment does not contain. That calculation has been removed.')
    print('  The experiment establishes that DPP6 is required for the complex to ACCUMULATE, and')
    print('  sits upstream of Kv4 and KChIP3 in assembly - consistent with DPP6 also being the')
    print('  longest-lived subunit of the mature complex. It says nothing about mature-state')
    print('  exchange, so it is out of the turnover argument.')
    print('  Note also DIV 10 here vs DIV 18-19 in Dorrbaum: different developmental stages.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    pd.DataFrame(rows).to_csv(f'{OUT}/tables/er_triage_requirement.csv', index=False)
    print(f'\nwrote {OUT}/tables/er_triage_requirement.csv')

    print('\n' + '=' * 100)
    print('WHERE THE EXCHANGE INFERENCE STANDS AFTER BOTH OBJECTIONS')
    print('=' * 100)
    print('  Still standing, because it never depended on the RNAi: within one preparation, KChIP4')
    print(f'  (2.02 d) and KChIP1 (2.67 d) turn over faster than Kv4.3 ({T_COMPLEX} d), and a mixture')
    print('  whose protected component decays at the slower rate cannot present the faster apparent')
    print('  rate. ER triage cannot rescue a protected-pool model, per objection 1.')
    print('  Withdrawn: any use of the lab\'s RNAi as evidence about mature-complex turnover. It is')
    print('  an assembly experiment, and it is now out of the turnover argument entirely.')
    print('  Unchanged: the co-IP test is the decisive one. What fraction of KChIP4 is bound to Kv4?')


if __name__ == '__main__':
    main()
