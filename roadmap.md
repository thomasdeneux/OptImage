
### OPTIMAGE BUGS

cleared list!!

### OPTIMAGE IMPROVE

**display**

* allow resize of full window
* put a label to image similar to time courses
* put real condition name rather than C0, C1… as much as possible
* Remove the small button and replace TRIALS/REGIONS info by additional info in the time courses title

**function**

* when opening files, read in chronological order rather than alphabetical
* forbid data/datop and signal/signalop at the same time, or really display 4 curves!!
* improve channel scale
* more spatial displays:
 - add trialcondition and conditiontrial
 - make 'update/autoupdate' also available for 'full options...', at least when it makes sense

**fn_4Dtoolbox**

* 'externalize' the edition of shapes, based on a class similar to interactivePolygon but even more generic; then a 'spline' mode can be easily added

**ergonomy**

* remove deprecated signal operations
* check for unsaved tpv file upon closing
* file panel: how files are sorted and how they are grouped
* Add a help menu that open the pdf of the powerpoint
* 'about' with the version of OptImage
* add a tutorial?
* build a displays manager to visualize/import/export program organization

**code organization**

* Make the function getsignals applicable not only to V but also to C: actually, make a big change by making signals getaccess private to content, and calling getsignals with the specifications (e.g. datop, signalop, all conditions)
* In tps_trial getcondition function, check whether ‘all conditions’ is stored already: in this case use it (maybe also do not compute/store the requested subcondition(s)?)

**packaging**

* make the code available (licensing in a separate process)
* consider making tpview-scripts an “outside” folder
* check that the default options have a nice set of displays


