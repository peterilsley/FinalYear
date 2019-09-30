A = gpuArray(rand(2^16,1));%transfers the number onto the GPU (SEND)

B = fft(A);
C=gather(B);%returns the data back into local memory (COME BACK)