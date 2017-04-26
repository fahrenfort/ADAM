function output = sample2time(samples,sample_rate)
    % samples to time
    output = fix((samples/sample_rate)*1000);