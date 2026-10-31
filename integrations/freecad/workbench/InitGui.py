"""Register HornSculpt's local geometry-link workbench."""

import FreeCADGui as Gui


class HornSculptWorkbench(Gui.Workbench):
    MenuText = "HornSculpt"
    ToolTip = "Geometrie mit HornSculpt verbinden"

    def Initialize(self):
        import hornsculpt_commands

        hornsculpt_commands.register_commands()
        commands = ["HornSculpt_Connect", "HornSculpt_Disconnect"]
        self.appendToolbar("HornSculpt", commands)
        self.appendMenu("HornSculpt", commands)

    def GetClassName(self):
        return "Gui::PythonWorkbench"


Gui.addWorkbench(HornSculptWorkbench())
