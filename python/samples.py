import numpy as np


def ranked_samples(N, compensate_slow_down=True):
    """gives a vector of length N, where the i'th entry corresponds to the
    expected race time of the runner who will start in the i'th position.
    Uses provided data as reprasentative of the fraction of people who will run
    a given time. Given data corresponds to 4 normal distribution (4 distinct
    groups of people) and tells us how many people will fall into each group.
      """
    np.random.seed(1)
    #PDF values [mu, sigma]
    wave1dist1 = [2015.52, 190.8]
    wave1dist2 = [2641, 417.5]
    wave2      = [3004.5, 596.5]
    wave3      = [3455.5, 835.5]
    if compensate_slow_down:
        wave2[0] = wave2[0] - 30
        wave3[0] = wave3[0] - 60
    wavefractions = [0.0715654, 0.150214, 0.494024, 0.284197]
    frac1      = wavefractions[0]
    frac2      = wavefractions[1]
    frac3      = wavefractions[2]
    frac4      = wavefractions[3]
    worldrecord = 1577.0
    distance   = 10000.0
    ulist = []

    N1 = int(np.floor(N*frac1))
    N2 = int(np.floor(N*frac2))
    N3 = int(np.floor(N*frac3))
    N4 = N - N1 - N2 - N3

    Ulist = [0.0]*N

    # Take N1 samples from distribution 1
    count = 0
    for i in range(0, N1):
        samp = 0
        while samp <worldrecord:
            #biased sampling, ignoring all sub record times, should be lognormal?
            samp = np.random.normal(wave1dist1[0], wave1dist1[1])
        Ulist[count] = distance/samp
        count +=1

    for i in range(0, N2):
        samp = 0
        while samp <worldrecord:
            #biased sampling, ignoring all sub record times, should be lognormal?
            samp = np.random.normal(wave1dist2[0], wave1dist2[1])
        Ulist[count] = distance/samp
        count +=1

    for i in range(0, N3):
        samp = 0
        while samp <worldrecord:
            #biased sampling, ignoring all sub record times, should be lognormal?
            samp = np.random.normal(wave2[0], wave2[1])
        Ulist[count] = distance/samp
        count +=1
    for i in range(0, N4):
        samp = 0
        while samp <worldrecord:
            #biased sampling, ignoring all sub record times, should be lognormal?
            samp = np.random.normal(wave3[0], wave3[1])
        Ulist[count] = distance/samp
        count +=1




    return Ulist
