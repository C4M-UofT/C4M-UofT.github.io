---
layout: page
title: C4M Schedule for the MD Program
hide: true
---

{% assign sorted_lectures = site.data.c4m-md | sort: 'date' %} 
{% include schedule_table.html lectures=sorted_lectures time='12-2:30 PM' %}