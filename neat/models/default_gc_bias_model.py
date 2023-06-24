"""
The following model is the default for GC bias analysis
"""

import numpy as np

# The fallowing is the model used by previous versions of NEAT
default_window_size = 50
"""
The index of the following array = the GC count in a 50 base window.

examples:
index 0 corresponds to all 50 base windows in the reference with 0 Gs or Cs
index 25 corresponds to all 50 base windows in the reference with exactly 25 Gs or Cs

The values are calculated bias of the coverage at each of those concentration levels.

examples:
default_gc_bias[25] is the average coverage for all 50 base windows with 25 Gs or Cs in the reference.

Results near 1 mean that this particular GC concentration showed no bias in the reference dataset.
"""
default_gc_bias = np.array([0.25948222342397126, 0.3962426640329114, 0.4887289131998323, 0.5849014221689215,
                            0.7279504143747476, 0.8314297610968944, 0.9027899565152305, 0.9430873717195684,
                            0.9647837065785595, 0.9804383369627002, 0.9867894312952692, 0.993331961245743,
                            0.9993103100274058, 0.998917223246418, 0.9993811476281264, 1.0012185663383528,
                            1.0108245084926835, 1.0061430923013108, 1.0482627430005476, 1.0553377046446224,
                            1.1361109625621573, 1.0035314993564972, 0.9911814369148244, 0.9860384899685335,
                            1.0310994134455291, 0.980876662808423, 0.9745982860372921, 0.968349841746828,
                            0.9615028768469494, 0.9525576655555664, 0.9408772817987261, 0.9282638187780413,
                            0.9087390607016832, 0.876892217753752, 0.8499714722795447, 0.8204915576791397,
                            0.7809959155032594, 0.7462713043523627, 0.7203840177228688, 0.6592433013681249,
                            0.6467600329348031, 0.5909430313687961, 0.5591719472549926, 0.5656766092749199,
                            0.469508521726588, 0.4173746102200429, 0.37560923897559, 0.2855614146018458,
                            0.2637100319267415, 0.2050777809758164, 1.0])
