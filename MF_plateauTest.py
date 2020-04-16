import numpy as np

def MovingAv(data: list,k=5):
    mask = np.ones(k)/k
    return np.convolve(data, mask, mode='same')

def rms_peak_ratio(data):
    rms = np.sqrt(np.mean(data**2))
    peak = np.max(np.abs(data))
    
    return rms, peak, rms/peak