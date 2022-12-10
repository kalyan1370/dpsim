---
title: "Concepts"
linkTitle: "Concepts"
weight: 4
---

The book introduces the reader to the general concepts implemented in DPsim, a dynamic phasor (DP) real-time simulator, as well as the physical models of the power system components that are used in simulations.
The first chapters give an overview of dynamic phasors and nodal analysis which are the two pillars of the main solver implemented in DPsim.
The second part describes in detail what are the physical equations for each model and how they are transformed and implemented for dynamic phasor simulations and other domains that are also supported by DPsim.

In order to be able to run a dynamic simulation, DPsim also includes a loadflow solver to compute the initial state of the network if it is not included in the network data.
Besides DP simulations, DPsim also comes with EMT models for some components which are used as reference for testing the DP models.
