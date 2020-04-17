import numpy as np

def MovingAv(data: list,k=5):
    mask = np.ones(k)/k
    return np.convolve(data, mask, mode='same')

def rms_peak_ratio(data):
    rms = np.sqrt(np.mean(data**2))
    peak = np.max(np.abs(data))
    
    return rms, peak, rms/peak

def fftspectrum(data, sampling_rate):
    time = np.linspace(start=1/sampling_rate, stop=len(data)*1/sampling_rate, num=len(data))
    freq = np.fft.fftfreq(time.shape[-1],d=1/sampling_rate)
    sp = np.fft.fft(data)
    
    return freq[(freq>0) & (freq<(sampling_rate/2))], sp[(freq>0) & (freq<(sampling_rate/2))]
    