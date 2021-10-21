# Tutorial 7: Temperature and Density Coupling of OpenMC, MOOSE, and THM

In this tutorial, you will learn how to:

- Couple OpenMC via temperature and density to separate MOOSE applications
  solving for the [!ac](T/H) physics in the solid and fluid domains
- Establish coupling between OpenMC and MOOSE for nested universe OpenMC models
- Apply homogenized temperature feedback to heterogeneous OpenMC cells

!alert! note
This tutorial makes use of the following major Cardinal classes:

- [OpenMCCellAverageProblem](/problems/OpenMCCellAverageProblem.md)

We recommend quickly reading this documentation before proceeding
with this tutorial. This tutorial also requires you to download a
mesh file from Box. Please download the files from the
`gas_assembly` folder [here](https://anl.app.box.com/folder/141527707499?s=irryqrx97n5vi4jmct1e3roqgmhzic89)
and place these files within the same directory structured
in `tutorials/gas_assembly`.
!alert-end!

In this tutorial, we couple OpenMC to the MOOSE heat conduction module
and [!ac](THM), a 1-D systems level thermal-fluids code based on
MOOSE. [!ac](THM) essentially contains all the single-phase physics in
RELAP-7 [!ac](relap7). OpenMC will receive temperature feedback from both
the MOOSE heat conduction module (for the solid regions) and from [!ac](THM)
(for the fluid regions). Density feedback will be provided by [!ac](THM).
This tutorial models a full-height [!ac](TRISO)-fueled prismatic gas reactor
fuel assembly. This application is
an extension of [Tutorial 6C](gas_compact.md), which coupled
OpenMC and MOOSE heat conduction (i.e. without an application providing fluid
feedback) for a unit cell version of a prismatic gas reactor assembly. In this
tutorial, we add fluid feedback and describe several nuances associated with
setting up feedback in OpenMC lattices.

This tutorial was developed with support from the NEAMS Thermal Fluids Center
of Excellence. A paper [!cite](novak2021)
describing the physics models and mesh refinement studies provides additional
context beyond the scope of this tutorial.

## Geometry and Computational Model

The geometry consists of a [!ac](TRISO)-fueled gas reactor assembly, loosely
based on a point design available in the literature [!ac](sterbentz).
A top-down view of the geometry is shown in [assembly].
The assembly is a graphite prismatic hexagonal block with 108 helium coolant channels,
210 fuel compacts, and 6 poison compacts. Each fuel compact contains [!ac](TRISO)
particles dispersed in a graphite matrix with a packing fraction of 15%.
All channels and compacts are ordered in a
triangular lattice with pitch $p_{cf}$.
Due to irradiation- and temperature-induced swelling of graphite, small helium gaps
exist between assemblies. In this tutorial, rather than explicitly model the inter-assembly
flow, we treat the gap regions as solid graphite.
There are also graphite reflectors above and below the assembly.

!media assembly.png
  id=assembly
  caption=[!ac](TRISO)-fueled gas reactor fuel assembly
  style=width:70%;margin-left:auto;margin-right:auto

The [!ac](TRISO) particles use a conventional design that consists of a central
fissil uranium oxycarbide kernel enclosed in a carbon buffer, an inner
[!ac](PyC) layer, a silicon carbide layer, and finally an outer
[!ac](PyC) layer. The geometric specifications for the assembly
dimensions are shown in [assembly], while dimensions for the
[!ac](TRISO) particles are summarized in [table1].

!table id=table1 caption=Geometric specifications for the [!ac](TRISO) particles
| Parameter | Value (cm) |
| :- | :- |
| TRISO kernel radius | 214.85e-4 |
| Buffer layer radius | 314.85e-4 |
| Inner PyC layer radius | 354.85e-4 |
| Silicon carbide layer radius | 389.85e-4 |
| Outer PyC layer radius | 429.85e-4 |

Heat is produced in the [!ac](TRISO) particles to yield a total power of
16.7 MWth. This heat is removed by helium flowing downwards through the coolant
channels with a total mass flowrate of 9.775 kg/s, which is assumed to be uniformly
distributed among the coolant channels. The outlet pressure is 7.1 MPa.

### Heat Conduction Model

### OpenMC Model

The OpenMC model is built using [!ac](CSG). The [!ac](TRISO) positions are sampled
using the [!ac](RSA) [algorithm in OpenMC](https://docs.openmc.org/en/stable/examples/triso.html).
OpenMC's Python [!ac](API) is
used to create the model with the script shown below. First, we define materials
for the various regions. Next, we create a single [!ac](TRISO) particle universe
consisting of the five layers of the particle and an infinite extent of graphite
filling all other space. We then pack pack uniform-radius spheres into a cylindrical
region representing a fuel compact, setting each sphere to be filled with the
[!ac](TRISO) universe.

Next, we loop over $n_l$ axial layers and create a unique hexagonal lattice for each layer.
This hexagonal lattice defines the fuel assembly structure, and consists of four different
universes: 1) a fuel pin plus surronding matrix (`f`), 2) a coolant channel plus surrounding matrix (`c`),
3) a boron carbide poision pin plus surrounding matrix (`p`), and 4) a homogeneous graphite hexagonal pincell to fill
the "boundaries" and centermost region (`g`). In each layer we set up the lattice
structure by listing the universes in each "ring" of the lattice, with `ring0` being
the centermost ring and `ring11` being the outermost ring. For example, `ring2` is the
innermost ring of fuel and coolant pins in this example, which is an alternation
of a fuel pin and a coolant pin.

Recall that temperatures in OpenMC can be set directly on the cell, but that fluid densities
can only be set on *materials*. For this reason, we need to create 108 unique coolant materials
for each axial plane if we want to be able to set unique densities in each coolant channel
region. Rather than creating 108 materials in a loop or through some other manual process,
we use the `clone()` feature in OpenMC to clone an existing coolant material 108 times per layer.
This duplicates the material properties (densities and isotopic compisition), but assigns
a new ID that allows individual tracking of density. The Python script used to create the
OpenMC model is shown below.



The level on which we will apply feedback from MOOSE is 1, because the geometry
consists of a hexagonal lattice (level 0), and we want to apply individual cell feedback
within that lattice (level 1). For the solid phase, this selection is equivalent to
applying a single temperature (per compact and per layer) for a compact region -
all [!ac](TRISO) particles and the surrounding matrix in each compact receives a uniform
temperature. Finally, to accelerate the particle tracking, we:

- Repeat the same [!ac](TRISO) universe in each axial layer and within each compact
- Superimpose a Cartesian search lattice in the fuel channel regions.

The OpenMC geometry, colored by either cell ID or instance, is shown in
[openmc_model]. Not shown are the axial reflectors on the top
and bottom of the assembly. The lateral faces are periodic, while the top
and bottom boundaries are vacuum.
The Cartesian search lattice in the fuel compact
regions is also visible in [openmc_model].

!media assembly_cells.png
  id=openmc_model
  caption=OpenMC model, colored by cell ID or instance
  style=width:60%;margin-left:auto;margin-right:auto

Cardinal applies uniform temperature and density feedback to OpenMC
for each unique cell ID $+$ instance combination. For this setup,
OpenMC receives on each axial plane a total of 721 temperatures and 108 densities
(one density per coolant channel). With references to the colors shown in
[openmc_model], the 721 cell temperatures correspond to:

\begin{equation*}
210\underbrace{\text{ fuel compacts}}_{\substack{\textup{1 TRISO compact (\textit{rainbow})}\\\textup{1 matrix region (\textit{purple})}}}\ \ +\ \ \ 108\underbrace{\text{ coolant channels}}_{\substack{\textup{1 coolant region (\textit{various})}\\\textup{1 matrix region (\textit{various})}}}\ \ +\ \ \ 6\underbrace{\text{ poison compacts}}_{\substack{\textup{1 poison region (\textit{brown})}\\\textup{1 matrix region (\textit{blue})}}}\ \ +\ \ \ 73\underbrace{\text{ graphite fillers}}_\text{1 matrix region (\textit{mustard})}
\end{equation*}

Because we will run OpenMC second, the initial fluid temperature is
set to a axial distribution given by bulk energy conservation ($q=\dot{m}C_{p,f}\left(T_f-T_{inlet}\right))
given the inlet temperature $T_{inlet}$, mass flowrate $\dot{m}$, fluid
isobaric specific heat $C_{p,f}$. Just for the purposes of obtaining a reasonable
fluid temperature initial condition, a sinusoidal heat source $q$ is assumed.
The fluid density is then set using the ideal gas [!ac](EOS) at a fixed pressure
of 7.1 MPa given the imposed temperature, i.e. $\rho_f(P,T_f)$.

To create the XML files required to run OpenMC, run the script:

```
$ python assembly.py
```

You can also use the XML files checked in to the `tutorials/gas_assembly` directory.

### THM Model

