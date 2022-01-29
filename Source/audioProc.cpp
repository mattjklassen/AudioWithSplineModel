/*
  ==============================================================================

    audioProc.cpp
    Created: 8 Jan 2022 3:57:29pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"


void MainContentComponent::changeState (TransportState newState)
{
    if (state != newState)
    {
        state = newState;

        switch (state)
        {
            case Stopped:
                playButton.setButtonText ("Play");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (false);
                transportSource.setPosition (0.0);
                break;

            case Starting:
                transportSource.start();
                break;

            case Playing:
                playButton.setButtonText ("Pause");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (true);
                break;

            case Pausing:
                transportSource.stop();
                break;

            case Paused:
                playButton.setButtonText ("Resume");
                stopButton.setButtonText ("Return to Zero");
                break;

            case Stopping:
                transportSource.stop();
                break;
        }
    }
}

void MainContentComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    if (playCycleOn) {
        currentSampleRate = sampleRate;
        updateAngleDelta();
    }
    if (playModelOn) {
        currentSampleRate = sampleRate;
    }
    if (playWavFileOn) {
        transportSource.prepareToPlay (samplesPerBlockExpected, sampleRate);
    }
}

void MainContentComponent::updateAngleDelta()
{
    auto cyclesPerSample = frequencySlider.getValue() / currentSampleRate;
    // cycles/sample = fractions of sampling frequency
    angleDelta = cyclesPerSample * 2.0 * juce::MathConstants<double>::pi;
    // angelDelta = frequency in radians/sample
}


    
void MainContentComponent::writeWavFile()
{
    if (sampleRate > 0)
    {
        // Create an OutputStream to write to our destination file...
        outputFile.deleteFile();
        if (auto fileStream = std::unique_ptr<FileOutputStream> (outputFile.createOutputStream()))
        {
            // Now create a WAV writer object that writes to our output stream...
            WavAudioFormat wavFormat;
            if (auto writer = wavFormat.createWriterFor (fileStream.get(), sampleRate, 1, 16, {}, 0))
            {
                fileStream.release(); // (passes responsibility for deleting the stream to the writer object that is now using it)
                if (writer != nullptr) {
                        writer->writeFromAudioSampleBuffer (writeBuffer, 0, writeBuffer.getNumSamples());
                }
                delete writer;
            }
        }
    }
}

void MainContentComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (playCycleOn) {
        auto* leftBuffer  = bufferToFill.buffer->getWritePointer (0, bufferToFill.startSample);
        auto* rightBuffer = bufferToFill.buffer->getWritePointer (1, bufferToFill.startSample);
        auto Pi = juce::MathConstants<double>::pi;
        auto localTargetFrequency = targetFrequency;
        auto frequencyIncrement = (localTargetFrequency - currentFrequency) / bufferToFill.numSamples;
//        int n = kVal + 3;
//        bcoeffRandomized += 1;
//        randBcoeff(bcoeffRandomized);
//        if (bcoeffRandomized == n) {
//            bcoeffRandomized = 0;
//        }
        for (auto sample = 0; sample < bufferToFill.numSamples; ++sample)
        {
            // compute currentSample value with spline of selected cycle:
            if (currentAngle > 2.0 * Pi) {
                currentAngle -= 2.0 * Pi;
//                control += 1;
//                if (control > 1) {
//                    control = 0;
//                }
            }
            control = 0;
            float currentSample = computeSpline(control, currentAngle / (2 * Pi));
            leftBuffer[sample]  = currentSample;
            rightBuffer[sample] = currentSample;
            sampleRendered += 1;
//            if (sampleRendered == samplesPerSelectedCycle - 1) {
//                sampleRendered = 0;
//            }

            currentFrequency += frequencyIncrement;
            currentAngle += angleDelta;
        }
    }
    if (playWavFileOn)
    {
        if (readerSource.get() == nullptr)
        {
            bufferToFill.clearActiveBufferRegion();
            return;
        }
        transportSource.getNextAudioBlock (bufferToFill);
    }
}


void MainContentComponent::readAudioData2 (AudioFormatReader *reader) {
    
    sampleCount = (int) reader->lengthInSamples;
    sampleRate = (int) reader->sampleRate;
    floatBuffer.setSize ((int) reader->numChannels, (int) reader->lengthInSamples);
    reader->read(&floatBuffer,                     // juce AudioBuffer <float>
                 0,                               // start sample in buffer
                 (int) reader->lengthInSamples,   // number of samples in file data
                 0,                               // start sample to fill in buffer
                 true,                            // use Left channel (0)
                 false);                          // use Right channel (1)
}

void MainContentComponent::readAudioData (File file) {
//    old version
//    juce::FileInputStream inputStream (file);
//    inputStream.read(buffer, 44);
//    dataSize = *reinterpret_cast<unsigned*>(buffer+40);
//    sampleRate = *reinterpret_cast<unsigned*>(buffer+24);
//    sampleCount = dataSize/2;
//    data = new char[dataSize];
//    inputStream.read(data, dataSize);
//    samples = reinterpret_cast<short*>(data);
//    std::cout << "size of wav file data is:  " << dataSize << std::endl;
//    std::cout << "sample rate from wav file is:  " << sampleRate << std::endl;
//    std::cout << "number of data samples is:  " << sampleCount << std::endl;
}


