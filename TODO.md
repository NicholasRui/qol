## Todo

Make some code for prop diagrams with different x axis options

make log spaced ticks that are good

Make hr diagram easily

- store crystallization curves

- store bagnulo & landstreet WDs

- store magnetic red giants

- store vrard, gaulme, yu catalogs

- store other special objects, eg GD 358, Zvrk, the fast-spinning RG from Tayar, GD 278

for mesatoolkit package:
- one possible step can be evolving a little bit to annihilate a certain element

- In Python package, artificially modify certain elements
make a custom mesa network that removes all hydrogen >:)

set directory for saving plots at top of doc and never touch it again

add easy capability for ticking by 2, 5, 10, or 1, 3, 10, etc.

add to mesatoolkit ability to change surface composition easily

add ability to write text parallel to plot

- include option to just read environment variables for MESA and mesasdk directories, without having to change config.py manually

TODO: Bug: Currently, if ./re fails, it does not terminate the script but simply goes to the next one, which may not need the one that failed. Then, if it saves a photo, it will never run the failed component. Make it so the products of a run MUST exist before moving on, or break.

TODO: add warning when there is no opacity setting specified -- since MESA will crash if this is not specified

Bug: when xlabel and ylabel are used before axis is initialized, it doesn't apply proper formatting to them

## automate animations?

incorporate kwargs idiomatically...
make a directory called "config" with text files for MESA edits, and move some of the code in config.py elsewhere

move some files like do_one, etc., into a dedicated folder
come up with a better directory name for pgplot


## automated workflow
make some "canonical" slides to mix and match
test points and stuff for figure development in package

In package, automatically plot magnetic white dwarfs B&L style

package should be able to import both histories and profile at the same time



# add something which automatically keeps track of fourier transforms for fft

- make a documentation website
- make this into a proper package
- requirements.txt

write a readme
- request people make PRs when they do things
- tell people what to cite
- nicholas rui



Draw plot locations for grid spec in ascii symbols


