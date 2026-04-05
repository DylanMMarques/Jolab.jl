````@raw html
---
layout: home

hero:
  name: Jolab
  text: Modelling framework for optical systems.
  tagline: Lorem ipsum...
  image:
    src: /logo.png
    alt: VitePress
  actions:
    - theme: brand
      text: Get Started
      link: /guide/what-is-vitepress
    - theme: alt
      text: View on GitHub
      link: https://github.com/vuejs/vitepress

features:
  - title: Modular modelling
    details: Divide complex optical systems into optical components. Simulate the complex system by propagating light through the optical components.
  - title: Library of optical solvers
    details: Each optical components provides its own optical solvers, enabling simulating of each component with the most appropriate solver, maximising the accuracy and speed.
  - title: High-performance
    details: Jolab uses KernelAbstractions.jl to be fast on both CPU and GPU (NVIDIA, AMD)
  - title: Differentiable
    details: Jolab is differentiable, enabling the use of gradient-based optimization algorithms for the design of optical systems.
---
````




