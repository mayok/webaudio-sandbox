function _range(start, stop, step = 1) {
  var index = 0;
  var length = (stop - start) / step;
  var result = Array(length);

  while(length) {
    result[index] = start;
    start += step;
    index += 1;
    length -= 1;
  }

  return result;
}

function mfcc( audioBuffer, { threshold = 0.01} ) {
  if (audioBuffer === null) return null;
  var sampleRate = 44100;
  var len = audioBuffer.length;

  // var center = len / 2;
  // var cuttime = 0.04;
  // var wav = audioBuffer.getChannelData(0);
  // var wavdata = wav.slice(center - cuttime/2*sampleRate, center + cuttime/2*sampleRate);

  var wavdata = preEmphHamming(audioBuffer.getChannelData(0));

  if( Math.max(...wavdata) < threshold )
    return null;

  // DFT
  // number of samples
  var n = 2048
  var im = [];
  [wavdata, im] = minifft(wavdata, n);


  var Pdft = [];
  for (var i = 0; i < n/2; i+=1) {
    Pdft[i] = wavdata[i] * wavdata[i] + im[i] * im[i]
  }

  // mel filterbank
  var numChannels = 20;
  var filterbank, fcenters;
  [filterbank, fcenters] = melFilterBank(sampleRate, n, numChannels);

  var mspec = [];
  _range(0, numChannels).map((c) => {
    mspec.push(
      Math.log10(
        filterbank[c].map((e, i) => Pdft[i] * e )
        .reduce((acc, val) => acc + val, 0)
      )
    );
  });

  return dct(mspec, norm=true).slice(0,12);
}

/*
* DCT: Discrete Cosine Transform
*
*
*            N-1
*  y[k] = 2* sum x[n]*cos(pi*k*(2n+1)/(2*N)), 0 <= k < N.
*            n=0
*/
function dct(signal, norm=false) {
  return _range(0, signal.length).map((k) => {
    return (norm ? (k===0 ? Math.sqrt(1 / (4 * signal.length)) : Math.sqrt(1 / (2 * signal.length))) : 1)
    * 2
    * (signal.map((x, n) => {
      return x * Math.cos( Math.PI * k * (2 * n + 1) / ( 2 * signal.length ))
    }).reduce((acc, val) => acc + val, 0));
  })
}

/*
* melFilterBank
*
*/
function melFilterBank(fs, nfft, numChannels) {
  fmax = fs / 2;
  melmax = 1127.01048 * Math.log(fmax / 700.0 + 1.0);
  nmax = nfft / 2;

  df = fs / nfft;
  dmel = melmax / (numChannels + 1)

  melcenters = _range(1, numChannels + 1).map((e) => e * dmel);
  fcenters = melcenters.map((e) => 700.0 * (Math.exp(e / 1127.01048) - 1.0));
  indexcenter = fcenters.map((e) => Math.round(e / df));

  indexstart = [0].concat(indexcenter.slice(0, numChannels - 1));
  indexstop  = indexcenter.slice(1, numChannels).concat([nmax]);

  filterbank = [];
  _range(0, numChannels).forEach((c) => {
    filterbank[c] = Array(nmax).fill(0);

    increment = 1.0 / (indexcenter[c] - indexstart[c]);
    _range(indexstart[c], indexcenter[c]).map((i) => {
      return filterbank[c][i] = (i - indexstart[c]) * increment;
    });

    decrement = 1.0 / (indexstop[c] - indexcenter[c]);
    _range(indexcenter[c], indexstop[c]).map((i) => {
      return filterbank[c][i] = 1.0 - ((i - indexcenter[c]) * decrement);
    });
  });

  return [filterbank, fcenters];

}

/*
* hamming window
*   w(n) = 0.54 - 0.46 cos( 2 PI n / (M - 1) ), 0 <= n <= M - 1
*/
function hamming(n) {
  var w = [];
  for(var i = 0; i < n; i += 1) {
    w[i] = 0.54 - 0.46 * Math.cos( 2.0 * Math.PI * i / (n - 1));
  }
  return w;
}

/*
* pre emphasis filter
*   y(n) = x(n) - p x(n - 1)
*/
function preEmphasis(signal, p = 0.97) {
  var y = [];
  y[0] = signal[0];
  for(var i = 1; i < signal.length; i+=1){
    y[i] = signal[i] - ( p * signal[i - 1]);
  }

  return y;
}

function preEmphHamming(signal) {
  var n = signal.length;
  var y = [];
  y[0] = signal[0] * 0.08;

  for(var i = 1; i < n; i++) {
    y[i] = signal[i] - (0.97 * signal[i-1]);
    y[i] *= 0.54 - 0.46 * Math.cos( 2.0 * Math.PI * i / (n - 1));
  }

  return y;
}

function minifft(data, n = 0) {
  var N = n || data.length;
  if(data.length < n) {
    var re = Array.prototype.slice.call(data).concat(new Array(n - data.length).fill(0));
  } else {
    var re = data.slice(0, n);
  }
  var im = new Array(n).fill(0);

  for (var i = 0; i < N; i+=1) {
    for(var j = 0, h = i, k = N; k >>= 1; h >>= 1) {
      j = (j << 1) | (h & 1);
    }
    if (j > i) {
      re[j] = [re[i], re[i] = re[j]][0];
      im[j] = [im[i], im[i] = im[j]][0];
    }
  }
  for(var hN = 1; hN * 2 <= N; hN *= 2) {
    for(var i = 0; i < N; i+= hN * 2) {
      for(var j = i; j < i + hN; j+=1) {
        var cos = Math.cos(Math.PI * (j - i) / hN),
        sin = Math.sin(Math.PI * (j - i) / hN);
        var tre =  re[j+hN] * cos + im[j+hN] * sin,
        tim = -re[j+hN] * sin + im[j+hN] * cos;
        re[j+hN] = re[j] - tre;
        im[j+hN] = im[j] - tim;
        re[j] += tre;
        im[j] += tim;
      }
    }
  }

  return [re, im];
}

function onError() {
  console.log("error happen");
}
