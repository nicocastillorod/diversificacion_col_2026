---
layout: page
title: Instructores
description: Instructores Taller Colevol 2026
---

## Instructores

{% assign instructors = site.staffers | where: 'role', 'Instructor' %}
{% for staffer in instructors %}
{{ staffer }}
{% endfor %}
