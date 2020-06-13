Idiot Comedy Steps Towards Fast Classical Shor's
===

You want to have the [_Curb Your Enthusiasm_ theme](https://www.youtube.com/watch?v=Ag1o3koTLWM) running throughout reading this one.

The quantum portion of Shor's is the period finding algorithm. I may have implemented a fast (with respect to expanding the amplitudes, not with the sieves which are of course faster for this specific instance) but classical period finding algorithm for `f(x) = (7^x) mod 15`. This would have strong comedy value.

Maybe not! I have to go scale it up, see if there's any catches in the general case. There are also parts that annoy me which may attack the validity of this thing.

The deterministic qubits are pretty confusing. I talk of them at length and implement fast Fourier transform comparable to QFT [here](https://github.com/howonlee/deterministic-qubit). Like, this is not going to make sense at all without reading that and it'll barely make sense with. I am a barely-make-sense sort of person.

I took the gate structure from [here](https://arxiv.org/pdf/quant-ph/0112176.pdf). They have quite a few optimizations, some of which I took.

We get to skip the quantum fourier transform entirely (but I implemented a deterministic thing as fast as a quantum fourier transform in the previous repo linked above, if you want to see one) because we don't have a no-copying theorem. So we just copy and find multiple satisfying sets, and find period with GCD. Strong comedy value.

Imagined Q&A
---

### OK, I didn't understand anything you said there.

Two multiple-input multiple-output logic gates of interest (actually, many more) are almost solely used in quantum computing but are perfectly compatible with classical computing. These are CNOT:

CNOT(a, b) = (a, a XOR b)

and CCNOT (Toffoli):

CCNOT(a, b, c) = (a, b, (a AND b) XOR c)

Of interest is the fact that these logic gates are involutive (they have an inverse, and are their own inverse). Toffoli gate [can be reduced to NAND](https://en.wikipedia.org/wiki/Toffoli_gate#Universality_and_Toffoli_gate), and importantly, modular exponentiation can be made out of these and single-qubit gates.

So if we add the ancilla bits like the quantum computing folks add ancilla qubits and take care of a lot of other stuff that the quantum computing folks have sort of taken care of for us, it is alleged but not proven that modular exponentiation is straightforwardly invertible. I do have an implementation for 7 ** x mod 15, though.

### OK, I didn't understand the answer to the previous question, either.

In order for modern mainstream public key cryptography that underlays, for example, your bank account existing and not being accessible to random people, to exist, some mathematics has to not exist.

That is, there has to not be a way to invert some certain functions called one-way functions, and specifically the function called modular exponentiation, which is said to be a one-way function.

The mathematics we want to not exist may, unfortunately, exist. I am not certain of this, but it may. I have an implementation for a small case.

### Aren't there real known quantum lower bounds?

BQP depends upon black box and quantum algorithms seem to me to still separate for black box functions. However, Shor's talks of a non-black box function, modular exponentiation. With modular exponentiation, we can look in the logical gate structure of the function.

### Should I panic?

Right now, no. Strong chance that this is just crankery and you can ignore this.

If I or someone gets a fast period-finding for modular addition out with this method (and not the trivial classical method), maybe. Prepare a big red button to push in case everything goes wahoonie-shaped.

Modular multiplication, start sweating.

Modular exponentiation, you should probably go panic and press the big red button. Things have gone wahoonie-shaped.

If the Grover's thing in L-infinity norm works out, there is no need to press the big red button because society as it exists currently will not survive the month. It probably won't work out.

### Why didn't you properly finish this thing before you go out and shove it on github?

This way, in the case that the whole thing goes through, you get like a fair bit of warning instead of none. Also, it increases the comedy value and I am optimizing on comedy value.

Also, I gotta do this on weekends, what with the job-having, so I was getting sick and tired of this and wanted to poke at the Grover's thing.

### You're a crank.

Maybe! I don't know! Let's see!

Remaining Questions For Me
----

Is it the backwards way or the forwards way? `7^x mod 15` is not actually complicated enough for it to matter in period finding, looks like.

Why didn't the output state matter? Was that a fluke?

Is this going to keep on working as we get to arbitrary modular addition, multiplication, exponentiation? The real workings of the operations are really unclear because of all the optimizations the Vandersypen et al group did, so I need to whack at it with really ordinary implementation of modular addition next, I think.

The statement I had in the last repo about the measurements with deterministic qubits being interpretation gets dicey when we actually want to take advantage of entanglement and interference. It's actually finding satisfying interpretations. This seemed very scary to me but is turning out to be less so. We get to do the satisfying quickly and in a superposition state, because of the invertible and involutive character of CNOT and CCNOT (Toffoli) gates: just apply them to invert, just like the Gaussian elimination for XORSAT but with a series of whole linear operators instead of elementary row operations. That last sentence there might be bullshit, but it hasn't proved to be bullshit with finding the period on our simple case.

Right now this finds period on `(7^x) mod 15`. Despite the questions, and unlike with actual quantum systems, however, I see no dire impediment towards getting this to find periods for whatever, it's just that more complicated modular exponentiation gets more complicated [really, really fast](https://arxiv.org/pdf/1207.0511.pdf). But still at a polyspace and polytime pace.

Mea culpa, I'm just one man. I will go and implement that arbitrary period finding, so the next news will probably be a full implementation of period finding for arbitrary values. I don't know if there will be unseen complications and I don't know if the speedup will stay, but it's worth a try.

I am also probably going to work on Grover's and the weirdness that L-infinity norm implies in Grover's after this one, probably simultaneously with arbitrary period finding. Really this one gives me many lances to tilt at the great windmill with, I could also see if composing OR out of Toffoli gates removes that difficulty with the OR in number partition with propositional representation of the dynamic program table.

Implementation comes first, proof later. If the whole thing works out proof is going to be a bit trivial anyways.

This one took a week's vacation from work, don't expect any of the proposed developments soon. I have toiled harder and longer for comedy value before, however. Again thanks to the YOSPOS folks.
