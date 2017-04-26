function output = time2sample(time, sample_rate)
    output = fix((time/1000)*sample_rate);