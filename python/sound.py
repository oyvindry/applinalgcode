import math, wave, sys, os, subprocess
from numpy import *

max_amplitude = 2**15-1 # iinfo('int16').max if numpy >= 1.0.3

def filter_impl(t, x, bd_mode):
    szx = shape(x)
    N = szx[0]
    n = int(prod(szx[1:]))
    y =  reshape(x, (N, n))
    tlen = len(t); N0 = int((tlen - 1)/2)
    w = 0
    
    if bd_mode.lower() == 'symm':
        w = concatenate([ y[N0:0:(-1)], y, y[(N-2):(N - N0 - 2):(-1)] ])
    elif bd_mode.lower() == 'per':
        w = concatenate([ y[(N - N0):], y, y[:N0]])
    elif bd_mode.lower() == 'none' or bd_mode.lower() == 'bd':
        w = concatenate( (zeros(N0, n), y, zeros(N0, n)) )
    for k in range(n):
        z = convolve(t, w[:, k])
        y[:, k] = z[(2*N0):(len(z)-2*N0)]
    x[:] = reshape(y, szx)
# End filter_impl

def audiowrite(filename, x, fs):
    """
    Writes the array data to the specified filename.
    The array data type can be arbitrary as it will be
    converted to numpy.int16 in this function.
    """
    ofile = wave.open(filename, 'w')
    ofile.setsampwidth(2)
    ofile.setframerate(fs)
    if x.ndim == 1:
        ofile.setnchannels(1)
    else:
        m, n = shape(x)
        ofile.setnchannels(n)
        x=x.flatten()
    x=max_amplitude*x
    x=x.astype(int16)
    x=x.astype(uint16)
    ofile.writeframesraw(x.astype(uint16).tostring())
    ofile.close()

def audioread(filename):
    """
    Read sound data in a file and return the data as an array
    with data type numpy.float, together with the sampling rate.
    Each sound channel will be a column in the array.
    """
    ifile = wave.open(filename)
    channels = ifile.getnchannels()
    fs = ifile.getframerate()
    frames = ifile.getnframes()
    x = ifile.readframes(frames)
    x = fromstring(x, dtype=uint16) 
    x=x.astype(int16)
    x=x.astype(float)/max_amplitude
    soundx = x
    if channels > 1:
        soundx = x.reshape( ( int(len(x)/channels), channels ) )
    return soundx, fs

def get_status_output(*args, **kwargs):
    p = subprocess.Popen(*args, **kwargs)
    stdout, stderr = p.communicate()
    return p.returncode, stdout, stderr
    
def play(x, fs, player=None):
    """
    Play a file with array data.  (The array is first written to file
    using the audiowrite function so the data type can be arbitrary.)  The
    player is chosen by the programs 'open' on Mac and 'start' on
    Windows. On Linux, try different open programs for various
    distributions. If keyword argument `player` is given, only this
    spesific command is run.
    """
    tmpfile = 'tmp.wav'
    audiowrite(tmpfile, x, fs)

    if player:
        msg = 'Unable to open sound file with %s' %player
        if sys.platform[:3] == 'win':
            status = os.system('%s %s' %(player, tmpfile))
        else:
            status, output, err = get_status_output([player, tmpfile])
            msg += '\nError message:\n%s' %output
        if status != 0:
            raise OSError(msg)
        return

    if sys.platform[:5] == 'linux':
        open_commands = ['gnome-open', 'kmfclient exec', 'exo-open', 'xdg-open', 'open']
        for cmd in open_commands:
            status, output, err = get_status_output([cmd, tmpfile])
            if status == 0:
                break
        if status != 0:
            raise OSError('Unable to open sound file, try to set player'\
                              ' keyword argument.')

    elif sys.platform == 'darwin':
        get_status_output(['open', tmpfile])
    else:
        # assume windows
        os.system('start %s' %tmpfile)