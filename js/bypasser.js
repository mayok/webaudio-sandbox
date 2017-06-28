/**
 * @class bypasser
 * @extends audioWorkletProcessor
 */
class Bypasser extends audioWorkletProcessor {
  constructor() {
    super();
  }

  process(input, output) {
    let inputChannelData = input.getChannelData(0);
    let outputChannelData = output.getChannelData(0);
    outputChannelData.set(inputChannelData);
  }
}

registerProcessor('bypasser', Bypasser);
