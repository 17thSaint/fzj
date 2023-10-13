import Pkg
allpacks = split(ARGS[1],",")
for pack in allpacks
	Pkg.add(pack)
	println("Added $pack")
end
