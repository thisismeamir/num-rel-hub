from manim import *
from wrapper import *
class ContentScene(Scene):
    global wrapper
    def construct(self):
                # Background Color

        config.frame_width = 16  # Set new width
        config.frame_height = 9  # Set new height
        self.camera.frame_width = 16
        self.camera.frame_height = 9

        height = config.frame_height
        width = config.frame_width
        self.add(Rectangle(width=16, height=9, color=BLACK, fill_opacity=1).move_to(ORIGIN))
        
        # "Overview" Heading
        overview_heading = Text("FOLIATION", font_size=160, font="IBM Plex Mono", color=WHITE)
        self.play(FadeIn(overview_heading), run_time=0.5)
        self.wait(1)
        self.play(overview_heading.animate.shift(LEFT), run_time=0.5)
        self.wait(1)
        text = wrapper("Einsteins\' theory of general relativity holds all the dimensions, be it time, or space dimensions, with the same regard, making any distinction between space and time coordinates more a matter of convention than a fundamental aspect of the theory. However, both our experiences and the laws of physics—especially on large scales—indicate that differentiating time from spatial coordinates is the most natural approach for describing physical processes[1].",54)
        mob=[]
        for t in text:
            mob.append(Text(t, fill_opacity=1, font_size=198,font="IBM Plex Sans", color=WHITE).move_to(ORIGIN, aligned_edge=LEFT))
        if len(mob)==1:
            self.play(FadeIn(mob[0]))
            self.wait(1)
        else:
            for i in range(0,len(mob)):
                self.play(
                    FadeIn(mob[i]),mob[i].animate.shift((len(mob) - 2 * i)/len(mob) *1.8*UP), run_time=0.75
                )
        self.wait(5)
