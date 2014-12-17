PROJECTS = projects/solver projects/mesh projects/prepare

COMMANDS = clean install strip 

.PHONY: projects $(PROJECTS) $(COMMANDS)

projects: $(PROJECTS)

$(PROJECTS):
	cd $@ && $(MAKE)

$(COMMANDS):
	for d in $(PROJECTS); do $(MAKE) --directory=$$d $@; done
