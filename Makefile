SRC_DIRS	:= src fosito particle_core
CREATE_DIRS	:= output

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(CLEANSRC)

all: 		$(SRC_DIRS) | $(CREATE_DIRS)

clean: 		$(CLEANSRC)

$(SRC_DIRS):
		$(MAKE) -C $@
$(CREATE_DIRS):
		mkdir -p $@
$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

# Dependecy information
src:		fosito particle_core | $(CREATE_DIRS)
particle_core: 	fosito

