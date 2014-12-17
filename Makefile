PROJECTS = projects/solver projects/mesh projects/prepare

.PHONY: projects $(PROJECTS)

projects: $(PROJECTS)

$(PROJECTS):
	cd $@ && $(MAKE)

