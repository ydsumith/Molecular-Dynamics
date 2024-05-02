#The py extension has nothing to do with. It is for formatting only.

#1) Below code removes type 3 atoms from the system
atomsk newlammpsdata.lmp -select prop type 3 -rmatom select thinmodel.lmp
